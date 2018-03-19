set output 'Lap.png'
set term png  
set view map
unset key
unset surface
 
#set style textbox opaque margins  0.5,  0.5 noborder

set cntrparam bspline
set contour base
#set cntrparam levels discrete 1e-90,1e-70,1e-60, 1e-58, 1e-56, 1e-55, 1e-54
set cntrparam levels auto

set pm3d 
set lmargin at screen 0.1
set rmargin at screen 0.8
set bmargin at screen 0.2
set tmargin at screen 0.9

set xlabel "Segund Pico [s]"
set ylabel "Amplitute"

splot   "Lap.dat"  w  l 
 