 gnuplot -persist << EOF
 set term wxt
 set title 'Figk.dat'
 set pm3d map
 set size square
 set xrange [-3.1416:3.0788]
 set yrange [-3.1416:3.0788]
 splot 'Figk.dat'
 #set term png size 1920,1280
 #set out 'Figk.dat.png'
 #rep
 #
EOF
 gnuplot -persist << EOF
 set term wxt
 set title 'Figk.dat'
 set nokey
 set grid
 set view 50,10,1,1
 splot 'Figk.dat' with pm3d
 #set term png size 1920,1280
 #set out 'Figk.dat.png'
 #rep
EOF
