# CBC tests

An extensive set of tests for the COIN-OR CBC MIP Solver. The main purpose is
to provide an easy to extend set of tests for the C Interface, and consequently
for CBC, that uses the github actions infrastructure.

Organization:

- instances/*.mps.gz: test instances
- instances/instances.csv:  a CSV file with the following data,
  instanceFile,relaxObj,bestBound,incumbentSolCost,optimal
