#! /bin/bash

ml python
grep "^Compression time" $1 | grep -o [0-9][0-9]*\\.[0-9]* > tmp
python3 compute_perf_all.py 
mv perf.txt comp_perf.txt
rm tmp
grep "^Decompression time" $1 | grep -o [0-9][0-9]*\\.[0-9]* > tmp
python3 compute_perf_all.py
mv perf.txt decomp_perf.txt
rm tmp

