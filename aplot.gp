set output 'alpha.png'
set term png  

set lmargin at screen 0.2
set rmargin at screen 0.8
set bmargin at screen 0.2
set tmargin at screen 0.9

set xrange [0.1:30] 
set xlabel "alpha"
set ylabel "P"
set logscale x
plot   "Lalpha.dat"  w  l notitle
unset logscale x

 
set output 'temp.png'
set term png  

set lmargin at screen 0.2
set rmargin at screen 0.8
set bmargin at screen 0.2
set tmargin at screen 0.9

set xlabel "Temperature"
set ylabel "P"
set xrange [2:8] 

plot   "Ltemp.dat"  w  l notitle


set output 'tau.png'
set term png  

set lmargin at screen 0.2
set rmargin at screen 0.8
set bmargin at screen 0.2
set tmargin at screen 0.9

set xrange [0.1:15] 
set xlabel "tau"
set ylabel "P"

plot   "LTau1.dat"  w  l notitle

set output 'tp.png'
set term png  

set lmargin at screen 0.2
set rmargin at screen 0.8
set bmargin at screen 0.2
set tmargin at screen 0.9

set xrange [0.1:15] 
set xlabel "T brust"
set ylabel "P"

plot   "Ltp.dat"  w  l notitle
