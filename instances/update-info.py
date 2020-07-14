from sys import argv, exit
from mip import Model, OptimizationStatus

ifile = argv[1]
iname = ifile.replace(".mps.gz", "")

objr = ""

insts = set()

f = open("instances.csv", "r")
for l in f:
    cols = l.split(",")
    inst = cols[0]
    insts.add(inst)
f.close()

if iname in insts:
    print("info for instance {} already computed.".format(iname))
    exit()

# relaxation info
m = Model(solver_name="gurobi")
m.read(ifile)
m.optimize(relax=True)
if m.status == OptimizationStatus.OPTIMAL:
    objr = m.objective_value
elif m.status == OptimizationStatus.INFEASIBLE:
    objr = "inf"
elif m.status == OptimizationStatus.UNBOUNDED:
    objr = "unb"
else:
    exit(1)

# mip info
m = Model(solver_name="gurobi")
m.read(ifile)
m.optimize(max_seconds=8000)
if m.status == OptimizationStatus.OPTIMAL:
    objmip = m.objective_value
    objbnd = m.objective_bound
elif m.status == OptimizationStatus.INFEASIBLE:
    objbnd = "inf"
    objmip = "inf"
else:
    exit(1)

f = open("instances.csv", "a")
f.write(
    "{},{},{},{},{}\n".format(
        iname, objr, objbnd, objmip, m.status == OptimizationStatus.OPTIMAL
    )
)
f.close()

if m.status in [OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE]:
    m.write("{}.mst".format(iname))
