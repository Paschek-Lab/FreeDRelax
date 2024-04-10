#!/bin/bash

#!/bin/bash


cat>diffu.inp<<EOF 
1  4.0 0.004 1.0 1.0 1000000
EOF
../../src/diffu < diffu.inp

cat norm_g2_of_t.dat | ../../src/scale_freed_hwang -sfac 2.53 -D 1.0 -d 1.0 -R 4.0 > norm_g2_cor_hf.dat
cat norm_g2_cor_hf.dat | awk '{print $1, $2-$3}' > diff.dat
