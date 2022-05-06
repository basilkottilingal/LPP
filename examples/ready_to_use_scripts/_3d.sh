#!/bin/bash
for u in particle-*;
do
	filename="${u%.*}"
	number=`echo $filename | sed 's/particle-//g';`
	echo "$number"
	gnuplot -persist <<-EOFMarker
	set title "fe forces" font ",14" textcolor rgbcolor "royalblue"
	set view 60,75
	set term png size 1024,1024
	set output "img-$number.png"
	set xlabel 'x'
	set ylabel 'z'
	set zlabel 'y'
	set xrange [-0.05:0.05]
	set yrange [-0.05:0.05]
	set zrange [-0.05:0.05]
	set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )
	#plot '$u' u 1:2 pt 7 lt 2 lc rgb 'red' t 'particles'
	splot '$u' u 1:3:2:4 pt 7 lt 2 lc palette z
	EOFMarker
done

