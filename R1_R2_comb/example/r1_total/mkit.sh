#!/bin/bash



paste ../r1_inter/R1_inter_298.dat ../r1_intra/R1_intra_298.dat | grep -v '#'  | awk '{print $1, $2+$10}'> R1_total_298.dat
paste ../r1_inter/R1_inter_273.dat ../r1_intra/R1_intra_273.dat | grep -v '#'  | awk '{print $1, $2+$10}'> R1_total_273.dat

paste ../r1_inter/R1_inter_298.dat ../r1_intra/R1_intra_298.dat | grep -v '#'  | awk '{print $1, $3+$11}'> R2_total_298.dat
paste ../r1_inter/R1_inter_273.dat ../r1_intra/R1_intra_273.dat | grep -v '#'  | awk '{print $1, $3+$11}'> R2_total_273.dat

paste ../r1_inter/R1_inter_298_select.dat ../r1_intra/R1_intra_298_select.dat | grep -v '#'  | awk '{print $1, $2+$10}'> R1_total_298_select.dat
paste ../r1_inter/R1_inter_273_select.dat ../r1_intra/R1_intra_273_select.dat | grep -v '#'  | awk '{print $1, $2+$10}'> R1_total_273_select.dat

paste ../r1_inter/R1_inter_298_select.dat ../r1_intra/R1_intra_298_select.dat | grep -v '#'  | awk '{print $1, $3+$11}'> R2_total_298_select.dat
paste ../r1_inter/R1_inter_273_select.dat ../r1_intra/R1_intra_273_select.dat | grep -v '#'  | awk '{print $1, $3+$11}'> R2_total_273_select.dat
