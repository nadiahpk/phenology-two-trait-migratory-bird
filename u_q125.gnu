# Trying to plot the mismatch etc as filled regions
name = "u_q125"
set terminal postscript eps size 5,5 enhanced color \
    font 'Helvetica,40' linewidth 2
set output name.'.eps'

set multiplot layout 2, 1
set lmargin at screen 0.2
set rmargin at screen 0.8

set xrange [135:175] # mod this

# Plot xc v n
#set yrange [0:300] # mod this
set ytics 2,1 # mod this
set xtics format ""
set xtics 10

set tmargin at screen 0.9
set bmargin at screen 0.7
set xtics
# set ylabel 'Population'
plot name.'.dat' u 1:4 w lines lw 4 lc rgb "blue" notitle

# Plot xc v arrival time etc
# set ylabel 'Key times'
set yrange [120:160] # mod this
set xtics format "%g"
set ytics 10
set tmargin at screen 0.7
set bmargin at screen 0.1
set style fill solid 
set xtics

plot name.'.dat' u 1:2:($2+$3) w filledcurves lc rgb "#DCDCDC" notitle ,\
     name.'.dat' u 1:($2+$3):1 w filledcurves lc rgb "#7F7F7F" notitle ,\
     name.'.dat' u 1:1:($2+$3+(14)) w filledcurves lc rgb "black" notitle ,\
     name.'.dat' u 1:2 w lines lw 4 lt 1 lc rgb "magenta" notitle,\
     name.'.dat' u 1:($2+$3) w lines lw 4 lt 1 lc rgb "red" notitle,\
     name.'.dat' u 1:($2+$3+(14)) w lines lw 4 lt 1 lc rgb "green" notitle


unset multiplot

