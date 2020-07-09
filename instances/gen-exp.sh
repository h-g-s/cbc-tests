#!/usr/bin/bash

rm -f exps.sh

for file in *.mps.gz;
do
    echo "pypy3 $file" >> exps.sh
done


