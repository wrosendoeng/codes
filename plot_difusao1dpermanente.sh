#! /bin/bash
#! Para compilar, basta escrever no WSL2: sed -i -e 's/\r$//' plot_difusao1dpermanente.sh
./difusao1dpermanente.exe
#! ./debug_difusao1dpermanente.exe

frame ( )
  {
    echo "set terminal pngcairo size 960,640 font 'Arial, 14' enhanced"
    echo "set termoption enhanced"
    echo "set encoding iso_8859_1"
    echo "set output 'grafico-psi5-alfa12-100pontos.png'"
    echo "set mxtics 5"
    echo "set mytics 5"
    echo "set grid mxtics mytics"
    echo "set xrange [0:1.0]"
    echo "set yrange [0:1.0]"
    echo "set xlabel 'X = x/L'"
    echo "set ylabel '{/Symbol q} = {/Symbol q}(X)'"
    echo "set key right box height 2"
    echo "set title 'Comparação para {/Symbol q} ({/Symbol a} = 1.2, {/Symbol y}^2 = 25, n = 100 pontos)' font ',20'"
    echo "plot 'trabalho01_psi5_alfa12_100pontos.txt' u 3:4 w p pt 7 title 'numérico', \\
    'trabalho01_psi5_alfa12_100pontos.txt' u 3:5 w lines dt 6 lc rgb 'black' lw 1.5 title 'analitico'"
  }

frame | gnuplot