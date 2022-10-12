set term pdfcairo color size 16cm,6.5cm font 'Times-New-Roman,14'
set out 'Nallusamy2006val-charge.pdf'
source="/mnt/c/Users/fetc9/Рабочий стол/РАН/Теплоаккумуляторы/1dim_u_PCM/Верификация/experimentalData.dat"#"exp.lnk"#
set multiplot
set lmargin 2
set rmargin 1
set xlabel 'Time, min' offset graph 0,graph 0.02
set ylabel 'PCM temperature, K' offset graph 0.02,graph 0
unset key
array dd[5] = [1,2,4,3,5]
array cc[5] = [1,6,7,2,4]
array ww[5] = [2,3,3,3,4]
array seta[4] = ["[2]","[3]","[4]","[5]"]
set size 0.44,1
#set xtics 0.2
set origin 0.05,0
set yrange [305:]
plot for [i=1:4] source ind "[PCM Temperature Nallusamy" u (column(1)/60):(column(i+1)) w lp ps .6 dt 1 lc cc[i] lw 1 ti "{/Times-New-Roman-Italic ".i."}",\
for [i=1:4] source ind "[OpenFOAM PCM Temperature Nallusamy" u (column(1)/60):(column(i+1)) w l dt dd[i] lc cc[i] lw 3 ti "{/Times-New-Roman-Italic ".(i+4)."}"
set origin 0.535,0
set ylabel 'Water temperature, ^oC'
set key r c rev Left vertical maxrows 4
#set xrange [60:140]
#set xtics 0.25
plot for [i=1:4] source ind "[HTF Temperature Nallusamy" u (column(1)/60):(column(i+1)) w lp ps .6 dt 1 lc cc[i] lw 1 ti "{/Times-New-Roman-Italic ".i."}",\
for [i=1:4] source ind "[OpenFOAM HTF Temperature Nallusamy" u (column(1)/60):(column(i+1)) w l dt dd[i] lc cc[i] lw 3 ti "{/Times-New-Roman-Italic ".(i+4)."}"
unset multiplot
set out

#------------------------------------------------------------------------------
