# plotD.plt
set encoding utf8
set terminal windows font "Consolas,12"
set nokey
set mouse
dados="Deslocamento.txt"
set zrange [:] reverse
set autoscale xfixmax
set autoscale xfixmin
set autoscale yfixmax
set autoscale yfixmin
set autoscale zfixmax
set autoscale zfixmin
set title "{/:Bold {/=14 Deslocamentos da Placa}}"
set xlabel "{/:Bold Eixo X [m]}" offset -5,-1
set ylabel "{/:Bold Eixo Y [m]}" offset 5,-1
set zlabel "{/:Bold Deslocamento [mm]}" rotate by 90 offset -2,0
set xtics 1
set mxtics 2
set ytics 1
set mytics 2
set grid xtics,mxtics,ytics,mytics,linetype rgb "#c0c0c0"
set contour base
set cntrparam bspline
set cntrlabel format start 10 interval -1 font 'Consolas,9'
set pm3d at s
set palette rgbformulae -10,-13,-26
splot dados using 1:2:3 with pm3d, dados using 1:2:3 with labels textcolor rgb "#46a2da", dados using 1:2:3 with lines linetype rgb "#77000000"