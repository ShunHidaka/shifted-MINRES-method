set terminal pdfcairo
set output "./Figure2.pdf"

PATH = "./stdout"
set grid
set xlabel "index of the shifted linear systems (m)"
set ylabel "number of iterations"
set yrange [0:]

# CLIQ6912std_A.csr
plot PATH."1_sminres.dat" u 1:4 w l lc 1 t "sMINRES",\
     PATH."1_scocg.dat"   u 1:4 w l lc 2 t "sCOCG"

# CLIQ55296std_A.csr
plot PATH."2_sminres.dat" u 1:4 w l lc 1 t "sMINRES",\
     PATH."2_scocg.dat"   u 1:4 w l lc 2 t "sCOCG"

# VCNT900h_A
plot PATH."3_sminres.dat" u 1:4 w l lc 1 t "sMINRES",\
     PATH."3_sbicg.dat"   u 1:4 w l lc 3 t "sBiCG"

# VCNT10800h_A
plot PATH."4_sminres.dat" u 1:4 w l lc 1 t "sMINRES",\
     PATH."4_sbicg.dat"   u 1:4 w l lc 3 t "sBiCG"
