"""Generates github workflow files to process tests."""

from random import shuffle
from glob import glob
from os import remove
from math import ceil

INSTS_PER_BATCH = 15

yml_header = """
name: C/C++ CI

on:
  push:
  pull_request:

jobs:
"""

# batch header
# x =
# """
#  test-batch-BATCHNUMBER:
#    runs-on: ubuntu-20.04
#    steps:
# """"

insts = []
ifeat = dict()

f = open("./instances/instances.csv")
for l in f:
    cols = l.split(",")
    iname = cols[0]
    ifeat[iname] = cols
    insts.append(iname)
f.close()

shuffle(insts)

wfiles = glob("./.github/workflows/*")
for wf in wfiles:
    remove(wf)

nWorkflows = ceil(len(insts) / INSTS_PER_BATCH)


f = open("./.github/workflows/test.yml", "w")
f.write(
    """
name: COIN-OR CBC Tests

on:
  push:
  pull_request:

jobs:
"""
)

for i in range(nWorkflows):
    j = 0
    f.write(
        """
  instances-batch-{}:
    runs-on: ubuntu-20.04
    steps:
      - uses: actions/checkout@v2
""".format(
            i
        )
    )
    while j < INSTS_PER_BATCH and insts:
        inst = insts.pop()
        c = ifeat[inst]
        f.write(
            """
      - name: test on {}
        run: ./bin/c-interface-solver ./instances/{}.mps.gz {} {} {} {}
""".format(
                inst, c[0], c[1], c[2], c[3], c[4]
            )
        )
        j += 1
f.close()
