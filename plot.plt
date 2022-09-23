unset logscale x
set xlabel "X" 
set ylabel "Y"
set grid
set yrange [-1:2]
set terminal png enhanced truecolor 
  set output 'tempfile.png'            
    plot "./dat1.txt" with lines title "Исходная функция", "./dat2.txt" with lines title "Стандартный многочленом", "./dat3.txt" with lines title "Многочлен Лагранжа"
  set out                              # restore the output redirection
set terminal GNUTERM                    