#! /bin/bash

dataset=$1
d_output=$1_"decomposition_time.txt"
touch $d_output
rm $d_output
grep Decomposition $dataset"_mgard.txt" | grep -o [0-9]*\\.[0-9]* >> $d_output
grep Decomposition $dataset"_opt_0.txt" | grep -o [0-9]*\\.[0-9]* >> $d_output
grep Decomposition $dataset"_opt_1.txt" | grep -o [0-9]*\\.[0-9]* >> $d_output
grep Decomposition $dataset"_opt_2.txt" | grep -o [0-9]*\\.[0-9]* >> $d_output
r_output=$1_"recomposition_time.txt"
touch $r_output
rm $r_output
grep Recomposition $dataset"_mgard.txt" | grep -o [0-9]*\\.[0-9]* >> $r_output
grep Recomposition $dataset"_opt_0.txt" | grep -o [0-9]*\\.[0-9]* >> $r_output
grep Recomposition $dataset"_opt_1.txt" | grep -o [0-9]*\\.[0-9]* >> $r_output
grep Recomposition $dataset"_opt_2.txt" | grep -o [0-9]*\\.[0-9]* >> $r_output

