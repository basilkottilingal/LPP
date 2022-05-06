list=system('ls -1B particle-*')
do for [fname in list] {
        cmd=sprintf("echo %s | sed 's/particle-//g'",fname)
        number=int(system(cmd))
        tname=sprintf("randomly moving particles, %d", number);
        set title tname font ",14" textcolor rgbcolor "royalblue"
        set term png size 1024,1024
        set output sprintf('img-%d.png',number)
        set xlabel 'x'
        set ylabel 'y'
        set xrange [0.:1.]
        set yrange [0.:1.]
        set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )
        plot fname u 1:2:3 pt 7 lt 2 lc palette z
}
