set output 'Lat.png'
set term png  
set view map
unset key
unset surface

 
#set style textbox opaque margins  0.5,  0.5 noborder

set cntrparam bspline
set contour base
#set cntrparam levels discrete 1e-90,1e-70,1e-60, 1e-58, 1e-56, 1e-55, 1e-54
set cntrparam levels discrete 0.25 ,0.5, 0.75

 
#set grid xtics nomxtics ytics nomytics noztics nomztics nortics nomrtics  nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
#set grid xtics nomxtics ytics nomytics noztics nomztics nortics nomrtics    nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set pm3d 
#set pm3d interpolate 0.1,1 flush begin noftriangles border linecolor rgb "gray"  linewidth 1.000 dashtype solid corners2color mean
 


set lmargin at screen 0.1
set rmargin at screen 0.8
set bmargin at screen 0.3
set tmargin at screen 0.9

  set yrange [2:8] 
  set xrange [0.2:15] 
  set logscale x
  #set xtics axis   0.1,10,15
 
  set xtics border out scale 1,1 mirror    offset character 0, -0.2, 0 autojustify
set ytics border out scale 1,1 

set xlabel "Alpha"
set ylabel "Temperatura"

 

splot   "Lat.dat"  w  l 
 