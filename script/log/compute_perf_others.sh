#! /bin/bash

ml python
grep "^Compression time" "$1"_"$2".log | grep -o [0-9][0-9]*\\.[0-9]* > tmp
python3 compute_perf.py $1
rm tmp

