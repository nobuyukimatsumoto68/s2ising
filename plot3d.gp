set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'
set output '3d.pdf'

# set xrange [-1.1:1.1]
# set yrange [-1.1:1.1]
# set zrange [-1.1:1.1]
set xrange [-0.8:1.1]
set yrange [-1.1:0.0]
set zrange [-0.8:1.1]
# set size ratio 1.0
set view equal xyz

set ytics
set xtics

set xlabel 'x'
set ylabel 'y'
set zlabel 'z'

# splot "simp1.dat", "dual1.dat", "duallink1.dat", "simplink1.dat", "mix1.dat"
splot "duallink1.dat", "simplink1.dat", "mix1.dat"