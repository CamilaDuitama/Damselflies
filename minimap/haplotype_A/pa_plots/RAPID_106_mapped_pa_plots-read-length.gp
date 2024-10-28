
        set ytics nomirror
        set y2tic nomirror
        set terminal png size 800,1000 truecolor
        set output "RAPID_106_mapped_pa_plots-read-length.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set grid noy2tics
        set multiplot layout 2,1 columnsfirst scale 1,0.9
        set yrange[0: 20]
        set ylabel "Read Count"
        set logscale x

        set xrange[50: 158]
        set xtics rotate by -45 out scale 1 nomirror
        set xlabel "Read Length"
        set boxwidth .8 relative
        set title "RAPID_106_mapped_pa.stats" noenhanced
        set title offset -15,-.5
        set key at graph 1, 1.15
        set y2range[0:2]
        plot "RAPID_106_mapped_pa_plots-read_len_rawdata.txt" using 1:2 lt rgb "orange" with boxes fs solid .7 title 'unbinned', "RAPID_106_mapped_pa_plots-read_len_histdata.txt" using 4:3 lt rgb "blue" with boxes title 'normalized bins(count/bin width)' axis x1y2
    
        set boxwidth 1 relative
        set y2range[0:20]
        plot "RAPID_106_mapped_pa_plots-read_len_rawdata.txt" using 1:2 lt rgb "orange" with boxes fs solid .7 title 'unbinned', "RAPID_106_mapped_pa_plots-read_len_histdata.txt" using 4:2 lt rgb "blue" with boxes title 'unnormalized bins' axis x1y2
    