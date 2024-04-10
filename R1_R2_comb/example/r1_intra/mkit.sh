#!/bin/bash

# Units
# -tau tau_KWW in ps
# -beta beta_KWW
# -A Pre factor A_KWW
# -rhh fixed H-H distance
# diff.dat: x: ps, y: 1
# diff.dat from system containing 8192 water molecules

#  rhh (TIP4P/2005) 0.1514
#  rhh (better)     0.1545

#rhh=0.1514
#rhh=0.1545
#rhh=1.6249
rhh=0.1561

../../src/gen_omega_MHz.pl >frequency_MHz.inp

head -501 diff_298.dat > diff.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.772 -tau 3.03 -beta 0.907 -dt 0.01 -ntmax 100000 -mhz -rhh $rhh  > Jomega_intra_298.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.772 -tau 3.03 -beta 0.907 -dt 0.01 -ntmax 100000 -mhz -relax -rhh $rhh  > R1_intra_298.dat

head -501 diff_273.dat > diff.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.799 -tau 6.30 -beta 0.884 -dt 0.01 -ntmax 100000 -mhz -rhh $rhh  > Jomega_intra_273.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.799 -tau 6.30 -beta 0.884 -dt 0.01 -ntmax 100000 -mhz -relax -rhh $rhh  > R1_intra_273.dat


cat >frequency_MHz.inp<<EOF
50.0
125.66
175.92
200.0
400.0
800.0
1200.0
EOF


head -501 diff_298.dat > diff.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.772 -tau 3.03 -beta 0.907 -dt 0.01 -ntmax 100000 -mhz -rhh $rhh  > Jomega_intra_298_select.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.772 -tau 3.03 -beta 0.907 -dt 0.01 -ntmax 100000 -mhz -relax -rhh $rhh  > R1_intra_298_select.dat

head -501 diff_273.dat > diff.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.799 -tau 6.30 -beta 0.884 -dt 0.01 -ntmax 100000 -mhz -rhh $rhh  > Jomega_intra_273_select.dat
cat frequency_MHz.inp | ../../src/sdens_kohlrausch -1H -A 0.799 -tau 6.30 -beta 0.884 -dt 0.01 -ntmax 100000 -mhz -relax -rhh $rhh  > R1_intra_273_select.dat


