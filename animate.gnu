set term pdfcairo lw 2 font "Times␣New␣Roman , 8 "
set output "4-order_RK_Method160_7.pdf"
set xlabel "t (s)"
set ylabel "z"
plot "4-order_RK_Method_more.dat" with lines title "4-order RK Method z versus time r=160.7"
set output
