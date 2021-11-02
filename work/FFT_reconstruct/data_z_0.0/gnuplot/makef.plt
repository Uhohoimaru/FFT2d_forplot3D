set terminal png
set xlabel "x"
set ylabel "f"
set key outside
set ticslevel 0

if (exist("n")==0 || n<0) n=0 
n=n+1
#set yr[-0.02:0.02]
strin = sprintf("../z_0_%05d",n)
strout = sprintf("./plot/all%03d.png",n)
gen = sprintf("m=%03d",n)
set output strout
#set logscale y
plot strin u 1:2 w l title gen


if (n<100)    reread

n=-1



system("convert -delay 10 -loop 0 ./plot/*.png ./reconstruct.gif")

