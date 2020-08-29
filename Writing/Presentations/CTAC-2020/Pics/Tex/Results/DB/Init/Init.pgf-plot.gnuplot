set table "Init.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set samples 2000; set dummy x; plot [x=-200:200] 2 -0.5*(1 + tanh(x/50));
