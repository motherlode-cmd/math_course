#     "a.txt" using 1:3 with lines t "mat"
#plot  "gauss_error.txt" u 2:3 with lines t "T(V) gauss analit"
#plot  "adams_error.txt" u 2:3 with lines t "V(t) adams_analit"
plot  "ralston_error.txt" u 2:3 with lines t "V(t) ralston analit"
#plot "a.txt" using 2:3 with lines t "x"

pause -1
