
            set terminal png size 600,500 truecolor
            set output "RAPID_106_mapped_pa_plots-gc-depth.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Mapped depth"
            set xlabel "Percentile of mapped sequence ordered by GC content"
            set x2label "GC Content [%]"
            set title "RAPID_106_mapped_pa.stats" noenhanced
            set x2tics ("30" 38.889,"40" 61.111,"50" 83.333)
            set xtics nomirror
            set xrange [0.1:99.9]

            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#dedede" t '10-90th percentile' , \
                 '-' using 1:2:3 with filledcurve lt 1 lc rgb "#bbdeff" t '25-75th percentile' , \
                 '-' using 1:2 with lines lc rgb "#0084ff" t 'Median'
        11.111	0.000	0.000
16.667	0.004	0.004
27.778	0.004	0.011
38.889	0.004	0.004
44.444	0.007	0.007
55.556	0.004	0.004
61.111	0.004	0.004
66.667	0.004	0.004
72.222	0.007	0.007
77.778	0.004	0.004
83.333	0.004	0.004
88.889	0.004	0.004
100.000	0.004	0.004
end
11.111	0.000	0.000
16.667	0.004	0.004
27.778	0.004	0.011
38.889	0.004	0.004
44.444	0.007	0.007
55.556	0.004	0.004
61.111	0.004	0.004
66.667	0.004	0.004
72.222	0.007	0.007
77.778	0.004	0.004
83.333	0.004	0.004
88.889	0.004	0.004
100.000	0.004	0.004
end
11.111	0.000
16.667	0.004
27.778	0.007
38.889	0.004
44.444	0.007
55.556	0.004
61.111	0.004
66.667	0.004
72.222	0.007
77.778	0.004
83.333	0.004
88.889	0.004
100.000	0.004
end
