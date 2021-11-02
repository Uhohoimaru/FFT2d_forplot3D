

if (exist("n")==0 || n<0) n=0 
n=n+1
strin_pm3d = sprintf("./FFT2d_gnuplot/data/FFT_2d_pm3d/fft2d_pm3d_z_%3.1f.txt",n)
strin_subave = sprintf("./FFT2d_gnuplot/data/subaverage/data_subaverage_z_%3.1f.txt",n)
strin_recon = sprintf("FFT2d_reconstruct_only_km_pm3d%3.1f.txt",n)


strout_pm3d_3d = sprintf("./FFT2d_gnuplot/gnuplot/all%03d.png",n)
strout_pm3d_2d = sprintf("./FFT2d_gnuplot/gnuplot/all%03d.png",n)
strout_subave  = sprintf("./FFT2d_gnuplot/gnuplot/all%03d.png",n)

gen = sprintf("m=%03d",n)

set output "./allplot.png"
set terminal pngcairo font "Times,16" enhanced size 1600,900
set multiplot layout 2,2

###############################################
set lmargin screen 0.1
set rmargin screen 0.4
set tmargin screen 0.9
set bmargin screen 0.6

set cblabel "PSD"
#set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
#set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
     0.8333 1 0 0, 1 1 1 1 )
#set term gif animate
set palette gray
 #set terminal postscript eps size 600,480 dashed enhanced color font 'Roman,10'


set pm3d 
set pm3d map
set xlabel "{k}"
set ylabel "{St_D}"
#set autoscale z
set logscale cb


#time =  sprintf('./FFT2dplot.png')
#set output time
#set title ttl
#set output time
set xr[-0.5:0.5]
set yr[0:8]
#set cbrange[100:200]
#set size ratio -1
 set pm3d interpolate 10,10
unset key
str = strin_pm3d
splot str with pm3d
unset logscale cb
#############################################################################################


unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 69, 333, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "k" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ -1.00000 : 1.00000 ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "StD" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ 0.00000 : 1.00000 ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ 0.00000 : 300.000 ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ 100.000 : 200.000 ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "ja_JP.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 10,10 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables
GNUTERM = "x11"
set lmargin screen 0.6
set rmargin screen 0.9
set tmargin screen 0.9
set bmargin screen 0.6

splot strin_pm3d with pm3d

############################################################################################



###########################################################################################
set lmargin screen 0.1
set rmargin screen 0.4
set tmargin screen 0.4
set bmargin screen 0.1

set cblabel ""
set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
#set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
     0.8333 1 0 0, 1 1 1 1 )
#set term gif animate
#set palette gray
 #set terminal postscript eps size 600,480 dashed enhanced color font 'Roman,10'
#set terminal pngcairo font "Times,14" enhanced size 800,450


set pm3d 
set pm3d map
set xlabel "{x}"
set ylabel "{t}"
set autoscale 


#time =  sprintf('./sub_ave.png')
#set output time
#set title ttl
#set output time
#set xr[0:40]
set yr[0:5]
#set cbrange[-0.001:0.001]
#set size ratio -1
# set pm3d interpolate 10,10
unset key
splot strin_subave with pm3d

#####################################################################################################



###############################################

set lmargin screen 0.6
set rmargin screen 0.9
set tmargin screen 0.4
set bmargin screen 0.1

set cblabel "p"
set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
#set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
     0.8333 1 0 0, 1 1 1 1 )
#set term gif animate
#set palette gray
 #set terminal postscript eps size 600,480 dashed enhanced color font 'Roman,10'


set pm3d 
set pm3d map
set xlabel "{x}"
set ylabel "{t}"
set autoscale 
#set logscale cb


#time =  sprintf('./FFT2dplot_reconstruct_kmins.png')
#set output time
#set title ttl
#set output time

#set xr[-0.5:0.5]
#set yr[0:8]
#set cbrange[0.1:1000]
#set size ratio -1
# set pm3d interpolate 10,10
unset key
#str = sprintf("./FFT2d_reconstruct_only_km_pm3d.txt")

splot strin_recon with pm3d

#############################################################################################

unset multiplot


if (n<31)    reread

n=-1
