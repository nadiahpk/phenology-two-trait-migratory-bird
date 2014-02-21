# Trying to plot the mismatch etc as filled regions
name = "u_q125"
set terminal postscript eps size 10,5 enhanced color \
    font 'Times,40' linewidth 8 dashlength 8
set output name.'.eps'

set style line 1 lt 2 lw 2 lc rgb "black" # pop size
set style line 2 lt 1 lw 3 lc rgb "blue" # hatching 
set style line 3 lt 5 lw 4 lc rgb "red" # laying date
set style line 4 lt 3 lw 4 lc rgb "brown" # arrival

set tics in
set tics scale 3,1

set multiplot layout 2, 2
set lmargin at screen 0.2
set rmargin at screen 0.5

set xrange [135:175] # mod this

# Plot xc v n
#set yrange [0:300] # mod this
set ytics 2,1 # mod this
set xtics format ""
set xtics 10

set tmargin at screen 0.975
set bmargin at screen 0.775
set xtics
set ylabel '{/Times-Italic n*/K}'
plot name.'.dat' u 1:4 w lines ls 1 notitle

# Plot xc v arrival time etc
# set ylabel 'Key times'
set yrange [120:160] # mod this
set xtics format "%g"
set ytics 10
set tmargin at screen 0.775
set bmargin at screen 0.175
set style fill solid 
set xtics

set xlabel '{/Times-Italic x_c}'
set ylabel 'Phenology' offset 1.3,0

plot name.'.dat' u 1:2:($2+$3) w filledcurves lc rgb "#DCDCDC" notitle,\
     name.'.dat' u 1:($2+$3):1 w filledcurve fill pattern 5 lc rgb "black" notitle ,\
     name.'.dat' u 1:1:($2+$3+(14)) w filledcurves lc rgb "#7F7F7F" notitle ,\
     name.'.dat' u 1:2 w lines ls 4 notitle,\
     name.'.dat' u 1:($2+$3) w lines ls 3 notitle,\
     name.'.dat' u 1:($2+$3+(14)) w lines ls 2 notitle

# legend -------------------------
t=0
set tmargin at screen 0.9
set bmargin at screen 0.1
set lmargin at screen 0.58
set rmargin at screen 0.95
set key center center
set noborder
unset tics
unset xlabel
unset ylabel
set yrange [0:1]
set xrange [0:1]
plot 2 t 'Population {/Times-Italic n*/K}' w lines ls 1, \
     2 t 'Hatching date {/Times-Italic x*}' w lines ls 2, \
     name.'.dat' u 1:($2+$3):1 w filledcurves lc rgb "#7F7F7F" lt 1 title 'Mismatch' ,\
     name.'.dat' u 1:2:($2+$3) w filledcurve fill pattern 5 noborder lt 1 lc rgb "black" title 'Incubation period' ,\
     2 t 'Laying date' w lines ls 3, \
     name.'.dat' u 1:1:($2+$3+(14)) w filledcurves lc rgb "#DCDCDC" lt 1 title 'Prelaying period {/Times-Italic z*}',\
     2 t 'Arrival date {/Times-Italic y*}' w lines ls 4


unset multiplot

