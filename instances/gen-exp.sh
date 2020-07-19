#!/usr/bin/bash

rm -f exps.sh

for file in *.mps.gz;
do
    echo "pypy3 update-info.py $file" >> exps.sh
done


