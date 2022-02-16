set ticslevel 0
set dgrid3d 30,30
set palette defined (0 "blue", 0.75 "white", 1.4 "red")
set style lines 100 lt 5 lw 0.5
set pm3d hidden3d 100
set grid
set view 74,216
unset key
splot 'fort.10' using 1:2:3 with pm3d
