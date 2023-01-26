### contour lines with labels
#reset session

f(x,y)=(x**2+y-11)**2+(x+y**2-7)**2

set xrange [0:5]
set yrange [0:5]
set isosample 250, 250

set contour base
set cntrparam levels disc 450,250,150,100,60,30,10,2
unset surface
set table $Contour
    splot f(x,y)
unset table

set style textbox opaque noborder

set multiplot layout 2,1
    plot $Contour u 1:2 w l lw 1.5 notitle, '' u 1:2:3 every 50 w labels boxed notitle

    set datafile commentschar " "
    plot for [i=1:8] $Contour u 1:2:(i) skip 5 index i-1 w l lw 1.5 lc var title columnheader(5)
unset multiplot
### end of code
