for i in *.dat ; do cp $i data ; gnuplot "plot1.plt" > $i.gif 2>/dev/null ; done
