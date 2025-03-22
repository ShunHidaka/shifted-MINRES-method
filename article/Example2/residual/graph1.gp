set terminal pdfcairo
set output "Figure1.pdf"

PATH = "./stderr"
set grid
set xlabel "number of iterations"
set ylabel "relative residual norm"
set logscale y
set format y "10^{%L}"

# CLIQ6912std_A.csr
plot PATH."1_sminres.dat" u 1:2 w l lc 1 t "sMINRES",\
     PATH."1_scocg.dat"   u 1:2 w l lc 2 t "sCOCG"

# CLIQ55296std_A.csr
plot PATH."2_sminres.dat" u 1:2 w l lc 1 t "sMINRES",\
     PATH."2_scocg.dat"   u 1:2 w l lc 2 t "sCOCG"

# VCNT900h_A
plot PATH."3_sminres.dat" u 1:2 w l lc 1 t "sMINRES",\
     PATH."3_sbicg.dat"   u 1:2 w l lc 3 t "sBiCG"

# VCNT10800h_A
plot PATH."4_sminres.dat" u 1:2 w l lc 1 t "sMINRES",\
     PATH."4_sbicg.dat"   u 1:2 w l lc 3 t "sBiCG"
