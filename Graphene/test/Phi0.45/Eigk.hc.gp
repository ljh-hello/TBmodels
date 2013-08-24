 gnuplot -persist << EOF
 set term wxt
 set title 'Eigk.hc'
 set pm3d map
 set size square
 set xrange [-3.1416:3.0788]
 set yrange [-3.1416:3.0788]
 splot 'Eigk.hc'
 #set term png size 1920,1280
 #set out 'Eigk.hc.png'
 #rep
 #
EOF
 gnuplot -persist << EOF
 set term wxt
 set title 'Eigk.hc'
 set nokey
 set grid
 set view 50,10,1,1
 splot 'Eigk.hc' with pm3d
 #set term png size 1920,1280
 #set out 'Eigk.hc.png'
 #rep
EOF
