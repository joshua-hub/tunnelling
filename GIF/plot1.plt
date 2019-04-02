!/usr/bin/gnuplot -persist

#  TITLES
set title "as12 melsom 42593249"
set xlabel "x Position"
set ylabel "Probability"

# AXES SETTINGS
set yrange [0:1]

set terminal gif
set output 'foobar.gif'

plot "" title= "B_0 = 50" w l 
