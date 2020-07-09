from mip import Model, OptimizationStatus
from glob import glob

insts = glob('*.mps.gz')

for inst in insts:
    iname = inst.replace('.mps.gz', '')

    objr = ''

    # relaxation info
    m = Model(solver_name='gurobi')
    m.optimize(relax=True)
    if m.status == OptimizationStatus.OPTIMAL:
        objr = m.objective_value
    elif m.status =
    

