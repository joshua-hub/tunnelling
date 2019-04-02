!/usr/bin/gnuplot -persist

#  TITLES
set title "as12 melsom 42593249"
set xlabel "x Position"
set ylabel "Probability"

# AXES SETTINGS
set yrange [0:1]
#set xrange [:1]
#set size square
#set xtics ('-2pi' -2*pi, '-pi' -pi, 0, 'pi' pi, '2pi' 2*pi) 
#set ytics 1
#set tics scale 0.75

# LEGEND
#set key default
#set key top left
#set key outside
#set key top left

# FILE OUTPUT
set term postscript color
set output "as11-problem3-melsom-42593249.ps"

f(x) = (50* exp(-16*(x**2)) + (x**2))/200

plot "temp50.dat" using 1:2 title "B_0 = 50" w l, \
     "temp100.dat" using 1:2 title "B_0 = 100" w l, \
     "temp150.dat" using 1:2 title "B_0 = 150" w l
