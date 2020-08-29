set table "Waves.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set samples 50; set dummy x; plot [x=0:11] 0.5 + 0.3*sin(pi*x/2) + 0.1*sin(pi*x/3) + 0.1*sin(pi*x/4) + 0.1*cos(pi*x/4) - 0.01*(x-5)**2 + 0.001*(x-5)**3 ;
