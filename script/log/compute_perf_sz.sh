#! /bin/bash

ml python
grep "^compression time" "$1"_sz.log | grep -o [0-9][0-9]*\\.[0-9]* > tmp
python3 compute_perf.py $1
rm tmp
