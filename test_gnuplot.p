# what i want to plot in the same code?
# downrange x height
# downrange x drift
# Mach x Cdrag
# Mach x Reynolds

# Main configuration:
set terminal jpeg
set output 'ballistics_trajectory.jpeg'
set datafile separator ','
set style data point
set style function lines
set size 1.0, 1.0
set multiplot
set grid
unset key

# Downrange x Height:
set size 0.5,0.5
set origin 0.0,0.5
set title "Downrange x Height"
set key below
plot \
    'out_nobb.csv' u 2:3 title "Inert base" with lines \
    'out_bb.csv' u 2:3 title "Base Bleed" with lines
set xlabel "Downrange (m)"
set ylabel "Height (m)" offset 1.0

# Downrange x Drift:
set size 0.5,0.5
set origin 0.0,0.0
set title "Downrange x Height"
set key below
plot \
    'out_nobb.csv' u 2:4 title "Inert base" with lines \
    'out_bb.csv' u 2:4 title "Base Bleed" with lines
set xlabel "Downrange (m)"
set ylabel "Drift (m)" offset 1.0

# Mach x CDrag:
set size 0.5,0.5
set origin 1.0,1.0
set title "Mach x Cdrag"
set key below
plot \
    'out_nobb.csv' u 9:8 title "Inert base" with lines \
    'out_bb.csv' u 10:9 title "Base Bleed" with lines
set xlabel "Downrange (m)"
set xrange [0:4]
set ylabel "Height (m)" offset 1.0
set yrange [0:0.5]

# Mach x CDrag
set size 0.5,0.5
set origin 1.0,0.5
set title "Mach x Reynolds"
set key below
plot \
    'out_nobb.csv' u 9:8 title "Inert base" with lines \
    'out_bb.csv' u 10:9 title "Base Bleed" with lines
set xlabel "Downrange (m)"
set xrange [0:4]
set ylabel "Height (m)" offset 1.0
set yrange [0:0.5]
