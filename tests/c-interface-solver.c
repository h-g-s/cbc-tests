/** C Interface test solver
 *
 * This test reads a MIP in .mps format, and performs the following tests:
 * - solve LP relaxation, objective value, solution feasibility
 * - solves MIP, checks best bound, solution cost, solution feasibility,
 *   optimality (if available)
 *
 * usage: c-interface-solver.c instanceFile relaxObj bestBound mipBound optimal
 *
 **/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Cbc_C_Interface.h>

/* instance name */
char instance_name[256] = "";

/* fill to inform current test step */
char test_step[512] = "";

/* fill to inform test element */
char test_element[512] = "";

#define ANNOUNCE_ERROR() { \
    fflush(stdout); fflush(stderr); \
    fprintf(stderr, "\n************************* ERROR *************************\n"); \
    fprintf(stderr, "  instance name: %s\n", instance_name); \
    fprintf(stderr, "      test step: %s\n", test_step); \
    fprintf(stderr, "   test element: %s\n", test_element); \
    fprintf(stderr, "    source code: %s:%d\n", __FILE__, __LINE__); \
    fflush(stdout); fflush(stderr); \
}

/* Checks for error if value is not exactly equal */
#define CHECK_DISCRETE_VALUE(VAR, VALUE) { \
    if (VAR != VALUE) { \
        ANNOUNCE_ERROR(); \
        fprintf(stderr, "        element: " #VAR "\n"); \
        fprintf(stderr, " expected value: %d\n", VALUE); \
        fprintf(stderr, "  current value: %d\n", VAR); \
        fprintf(stderr, "*********************************************************\n"); \
        fflush(stdout); fflush(stderr); success = 0; goto END; \
    }\
}

/* tolerances: larger than solver default tolerances,
 * focusing on more significant errors, which are more important
 * and easier to debug */
const double ABS_TOL = 1e-4;
const double REL_TOL = 0.01;

#define MIN(a, b) ((a<b) ? (a) : (b))
#define MAX(a, b) ((a>b) ? (a) : (b))

/* Checks for error if value is different with a tolerance, including relative tol  */
#define CHECK_BOUND_AT_LEAST(VAR, VALUE) {\
    {    double absdiff = MAX(ABS_TOL, fabs(VALUE * REL_TOL)); \
         if (VAR >= VALUE - absdiff ) { \
            ANNOUNCE_ERROR() ; \
            fprintf(stderr, "        element: " #VAR "\n"); \
            fprintf(stderr, " expected value: %g\n", VALUE); \
            fprintf(stderr, "  current value: %g\n", VAR); \
            fprintf(stderr, "     difference: %g\n", fabs(VAR)-VALUE); \
            fprintf(stderr, "*********************************************************\n\n"); \
            fflush(stdout); fflush(stderr); success = 0; goto END; \
        } } \
}

/* Checks for error if value is different with a tolerance, including relative tol  */
#define CHECK_BOUND_AT_MOST(VAR, VALUE) {\
    {    double absdiff = MAX(ABS_TOL, fabs(VALUE * REL_TOL)); \
         if (VAR <= VALUE + absdiff ) { \
            ANNOUNCE_ERROR() ; \
            fprintf(stderr, "        element: " #VAR "\n"); \
            fprintf(stderr, " expected value: %g\n", VALUE); \
            fprintf(stderr, "  current value: %g\n", VAR); \
            fprintf(stderr, "     difference: %g\n", fabs(VAR)-VALUE); \
            fprintf(stderr, "*********************************************************\n\n"); \
            fflush(stdout); fflush(stderr); success = 0; goto END; \
        } } \
}


/* Checks for error if value is different with a tolerance, including relative tol  */
#define CHECK_BOUND_EQUAL(VAR, VALUE) {\
    {    double absdiff = MAX(ABS_TOL, fabs(VALUE * REL_TOL)); \
         if (VAR <= VALUE-absdiff || VAR >= VALUE + absdiff ) { \
            ANNOUNCE_ERROR() ; \
            fprintf(stderr, "        element: " #VAR "\n"); \
            fprintf(stderr, " expected value: %g\n", VALUE); \
            fprintf(stderr, "  current value: %g\n", VAR); \
            fprintf(stderr, "     difference: %g\n", fabs(VAR)-VALUE); \
            fprintf(stderr, "*********************************************************\n\n"); \
            fflush(stdout); fflush(stderr); success = 0; goto END; \
        } } \
}

/* Checks for error if value is different with a tolerance  */
#define CHECK_CONTINUOUS_VALUE_EQUAL(VAR, VALUE) {\
    {   if (VAR <= VALUE-ABS_TOL || VAR >= VALUE + ABS_TOL ) { \
            ANNOUNCE_ERROR() ; \
            fprintf(stderr, "        element: " #VAR "\n"); \
            fprintf(stderr, " expected value: %g\n", VALUE); \
            fprintf(stderr, "  current value: %g\n", VAR); \
            fprintf(stderr, "     difference: %g\n", fabs(VAR)-VALUE); \
            fprintf(stderr, "*********************************************************\n\n"); \
            fflush(stdout); fflush(stderr); success = 0; goto END; \
        } } \
}



/* Checks for error if value is different with a tolerance  */
#define CHECK_CONTINUOUS_AT_MOST(VAR, VALUE) {\
    {   if (VAR >= VALUE + ABS_TOL ) { \
            ANNOUNCE_ERROR() ; \
            fprintf(stderr, "        element: " #VAR "\n"); \
            fprintf(stderr, " expected value: <= %g\n", VALUE); \
            fprintf(stderr, "  current value: %g\n", VAR); \
            fprintf(stderr, "     difference: %g\n", fabs(VAR)-VALUE); \
            fprintf(stderr, "*********************************************************\n\n"); \
            fflush(stdout); fflush(stderr); success = 0; goto END; \
        } } \
}

/* Checks for error if value is different with a tolerance  */
#define CHECK_CONTINUOUS_AT_LEAST(VAR, VALUE) {\
    {   if (VAR <= VALUE -ABS_TOL ) { \
            ANNOUNCE_ERROR() ; \
            fprintf(stderr, "        element: " #VAR "\n"); \
            fprintf(stderr, " expected value: >= %g\n", VALUE); \
            fprintf(stderr, "  current value: %g\n", VAR); \
            fprintf(stderr, "     difference: %g\n", fabs(VAR)-VALUE); \
            fprintf(stderr, "*********************************************************\n\n"); \
            fflush(stdout); fflush(stderr); success = 0; goto END; \
        } } \
}



int check_solution(Cbc_Model *m, char integrality, const double *x);

int main( int argc, char **argv ) {
    int success = 1;

    if (argc < 6) {
        fprintf(stderr, "usage: c-interface-solver instanceFile relaxObj bestBound mipBound optimal");
        exit(1);
    }

    printf("Starting test with the CBC C Interface. CBC Build info:\n\n");
    char binfo[2048];
    Cbc_getBuildInfo(binfo);
    printf("%s\n", binfo);

    char problemName[256] = "";
    strcpy(problemName, argv[1]);

    char *s = strstr(problemName, ".mps.gz");
    if (s)
        *s = '\0';

    s = (strstr(problemName, "instances/"));
    if (s)
        strcpy(instance_name, s+10);
    else
        strcpy(instance_name, s);

    double relaxObj = DBL_MAX;
    double bestBound = DBL_MAX;
    double mipBound = DBL_MAX;
    char opt = (strcasecmp(argv[5], "True") == 0);
    char relaxInf = (strcasecmp(argv[2], "inf") == 0);
    char mipInf = (strcasecmp(argv[4], "inf") == 0);
    if (!relaxInf) {
        relaxObj = atof(argv[2]);
    }
    if (!mipInf) {
        bestBound = atof(argv[3]);
        mipBound = atof(argv[4]);
    }

    Cbc_Model *m = Cbc_newModel();

    strcpy(test_step, "reading instance");
    Cbc_readMps(m, argv[1]);

    strcpy(test_step, "solving linear programming relaxation");
    Cbc_solveLinearProgram(m);

    strcpy(test_step, "checking solution of LP relaxation");
    if (relaxInf) {
        CHECK_DISCRETE_VALUE( Cbc_isProvenInfeasible(m), 1 );
        CHECK_DISCRETE_VALUE( Cbc_isProvenOptimal(m), 0 );
    }
    else {
        CHECK_DISCRETE_VALUE( Cbc_isProvenInfeasible(m), 0 );
        CHECK_DISCRETE_VALUE( Cbc_isProvenOptimal(m), 1 );
        CHECK_BOUND_EQUAL( Cbc_getObjValue(m), relaxObj );
    }

    if (!check_solution(m, 0, Cbc_getColSolution(m))) {
        success = 0;
        goto END;
    }

    Cbc_setIntParam(m, INT_PARAM_MAX_NODES, 1000);
    Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 300);

    strcpy(test_step, "integer optimization");
    Cbc_solve(m);

    strcpy(test_step, "checking optimization results");
    if (Cbc_isProvenInfeasible(m)) {
        CHECK_DISCRETE_VALUE(mipInf, 1);
    } else {
        if (Cbc_numberSavedSolutions(m)) {
            if (!check_solution(m, 1, Cbc_getColSolution(m))) {
                success = 0;
                goto END;
            }
            strcpy(test_step, "checking optimization status");
            CHECK_DISCRETE_VALUE(Cbc_isProvenInfeasible(m), 0);
            CHECK_DISCRETE_VALUE(Cbc_isAbandoned(m), 0);
            if (Cbc_isProvenOptimal(m)) {
                if (opt) {  // optimal cost available
                    strcpy(test_step, "checking bounds");
                    CHECK_BOUND_EQUAL(Cbc_getObjValue(m), mipBound);
                } else { // only bounds available
                    if (Cbc_getObjSense(m) == 1.0) { // minimize
                        CHECK_CONTINUOUS_AT_LEAST(Cbc_getObjValue(m), bestBound);
                    } else { // maximize
                        CHECK_CONTINUOUS_AT_MOST(Cbc_getObjValue(m), bestBound);
                    }
                } // checking bounds
            } // optimal
        } // solution found
    } // not infeasiblr

END:
    Cbc_deleteModel(m);
    if (success)
        return 0;

    return 1;
}

int check_solution(Cbc_Model *m, char integrality, const double *x) {
    char success = 1;

    char cname[256], str[256], rname[256];
    int j, nz, k;
    int numRows = Cbc_getNumRows(m);
    double lhs;
    const int *idx;
    const double *coef;
    double rhs, rlb, rub, obj_val = 0.0;

    if (integrality) {
        strcpy(test_step, "checking integrality of variables");
        for ( j=0 ; (j<Cbc_getNumCols(m)) ; ++j ) {
            if (Cbc_isInteger(m, j)) {
                Cbc_getColName(m, j, cname, 256);
                sprintf(test_element, "variable %s (%d)", cname, j);
                CHECK_CONTINUOUS_VALUE_EQUAL(x[j], floor(x[j]+0.5))
            }
            obj_val += Cbc_getColObj(m, j) * x[j];
        }
    } else {
        for ( j=0 ; (j<Cbc_getNumCols(m)) ; ++j ) {
            obj_val += Cbc_getColObj(m, j) * x[j];
        }
    }
    strcpy(test_element, "");

    strcpy(test_step, "checking computed objective value");
    CHECK_BOUND_EQUAL(Cbc_getObjValue(m), obj_val);

    strcpy(test_step, "testing if solution satisfies all problem constraints");
    for ( int i=0 ; (i<numRows) ; ++i ) {
        Cbc_getRowName(m, i, rname, 256);
        sprintf(test_element, "constraint %s (%d)", rname, i);
        lhs = 0.0;
        nz = Cbc_getRowNz(m, i);
        idx = Cbc_getRowIndices(m, i);
        coef = Cbc_getRowCoeffs(m, i);

        for ( k=0 ; (k<nz) ; ++k ) {
            lhs += x[idx[k]]*coef[k];
        }

        rhs = Cbc_getRowRHS(m, i);
        rlb = Cbc_getRowLB(m, i);
        rub = Cbc_getRowUB(m, i);

        switch (Cbc_getRowSense(m, i)) {
            case 'L':
                CHECK_CONTINUOUS_AT_MOST(lhs, rhs);
                break;
            case 'G':
                CHECK_CONTINUOUS_AT_LEAST(lhs, rhs);
                break;
            case 'E':
                CHECK_CONTINUOUS_VALUE_EQUAL(lhs, rhs);
                break;
            case 'R':
                CHECK_CONTINUOUS_AT_MOST(lhs, rub);
                CHECK_CONTINUOUS_AT_LEAST(lhs, rlb);
                break;
        }
    } // all rows
    strcpy(test_element, "");

    /* TODO : Test satisfaction of SOS */
END:

    return success;
}

