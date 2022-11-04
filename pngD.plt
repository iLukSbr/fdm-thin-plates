# pngD.plt
set encoding utf8
set terminal png font "Consolas,100" size 7680,4320 linewidth 2
set output "graf-desloc.png"
set nokey
dados="Deslocamento.txt"
set zrange [:] reverse
set autoscale xfixmax
set autoscale xfixmin
set autoscale yfixmax
set autoscale yfixmin
set autoscale zfixmax
set autoscale zfixmin
set title "{/:Bold {/=120 Deslocamentos da Placa}}"
set xlabel "{/:Bold Eixo X [m]}" rotate by -8 offset 0,-1
set ylabel "{/:Bold Eixo Y [m]}" rotate by 12 offset 3,-1
set zlabel "{/:Bold Deslocamento [mm]}" rotate by 90 offset -1,0
set xtics 1
set mxtics 2
set ytics 1
set mytics 2
set grid xtics,mxtics,ytics,mytics,linetype rgb "#c0c0c0"
set contour base
set cntrparam bspline
set cntrlabel format start 10 interval -1 font 'Consolas,85'
set pm3d at s
set palette rgbformulae -10,-13,-26
splot dados using 1:2:3 with labels textcolor rgb "#46a2da", dados using 1:2:3 with pm3d, dados using 1:2:3 with lines linetype rgb "#77000000"