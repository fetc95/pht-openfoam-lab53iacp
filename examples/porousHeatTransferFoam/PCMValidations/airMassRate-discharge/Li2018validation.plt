set term pdfcairo color size 16cm,6.5cm font 'Times-New-Roman,14'
set out 'Li2018val-PCM.pdf'
source="/mnt/c/Users/fetc9/Рабочий стол/РАН/Теплоаккумуляторы/1dim_u_PCM/Верификация/experimentalData.dat"#"exp.lnk"#
set multiplot
set lmargin 2
set rmargin 1
set xlabel 'Time, min'
set ylabel 'PCM temperature, K' offset graph 0.02,graph 0
unset key
array dd[5] = [1,2,4,3,5]
array cc[5] = [1,6,7,2,4]
array ww[5] = [2,3,3,3,4]
array seta[4] = ["[2]","[3]","[4]","[5]"]
set size 0.44,1
#set xtics 0.2
set origin 0.06,0
set xrange [0:60]
plot for [i=1:3] source ind "[PCM Temperature Li" u (column(1)/60):(column(i+2)) w lp ps 1 dt 1 lc cc[i] lw 1 ti "{/Times-New-Roman-Italic ".i."}",\
for [i=1:3] source ind "[OpenFOAM PCM Temperature Li" u (column(1)/60):(column(i+2)) w l dt dd[i] lc cc[i] lw 3 ti "{/Times-New-Roman-Italic ".(i+3)."}"
set origin 0.53,0
unset ylabel
set key r t rev Left vertical maxrows 3
set xrange [60:140]
#set xtics 0.25
replot
unset multiplot
set out

#------------------------------------------------------------------------------
