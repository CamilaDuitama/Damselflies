
            set terminal png size 600,400 truecolor
            set output "RAPID_106_mapped_pa_plots-gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "RAPID_106_mapped_pa.stats" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",29.65) at 29.65,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        11	0.000000
23	0.250000
25	0.000000
28	0.250000
28	0.500000
29	1.000000
30	0.750000
31	0.000000
33	0.250000
33	0.000000
34	0.250000
35	1.000000
36	0.750000
37	0.250000
38	0.000000
40	0.250000
44	0.000000
48	0.250000
49	0.500000
50	0.250000
51	0.000000
52	0.250000
54	0.000000
56	0.500000
58	0.000000
59	0.250000
end
