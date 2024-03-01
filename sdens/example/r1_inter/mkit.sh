#!/bin/bash
# Units
# -D inter diffusion coefficient: nm/ps
# -d distance of closest approach: nm
# -rho spin density: 1/nm**3
# -scale no units   (scale diff to infinite system size)
# diff.dat: x: ps, y: 1
# diff.dat from system containing 8192 water molecules


../../src/gen_omega_MHz.pl >frequency_MHz.inp

/bin/cp diff_273.dat diff.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H  -D 2.22e-3 -d 0.192 -mhz -scale 1.344 -rho 66.8367 > Jomega_inter_273.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H  -D 2.22e-3 -d 0.192 -mhz -relax -scale 1.344 -rho 66.8367 > R1_inter_273.dat

/bin/cp diff_298.dat diff.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H  -D 4.62e-3 -d 0.193 -mhz -scale 1.1967 -rho 66.663 > Jomega_inter_298.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H  -D 4.62e-3 -d 0.193 -mhz -relax -scale 1.1967 -rho 66.663 > R1_inter_298.dat


cat >frequency_MHz.inp<<EOF
50.0
200.0
400.0
800.0
1200.0
EOF

/bin/cp diff_273.dat diff.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H -D 2.22e-3 -d 0.192 -mhz -scale 1.344 -rho 66.8367 > Jomega_inter_273_select.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H -D 2.22e-3 -d 0.192 -mhz -relax -scale 1.344 -rho 66.8367 > R1_inter_273_select.dat

/bin/cp diff_298.dat diff.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H -D 4.62e-3 -d 0.193 -mhz -scale 1.1967 -rho 66.663 > Jomega_inter_298_select.dat
cat frequency_MHz.inp | ../../src/sdens_freed_hwang -1H -D 4.62e-3 -d 0.193 -mhz -relax -scale 1.1967 -rho 66.663 > R1_inter_298_select.dat

