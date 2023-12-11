#!/bin/bash


cat>diffu_rm2.inp<<EOF 
 0.001 100.0 1.0 1.0 1000000
EOF

../../src/diffu_rm2 < diffu_rm2.inp
cat norm_g2_of_t.dat | ../../src/scale_freed_hwang -D 1.0 -d 1.0 -R 100.0 > norm_g2_cor_hf.dat
cat norm_g2_cor_hf.dat | awk '{print $1, $2-$3}' > diff.dat
