#!/bin/bash

for i in MAplot_res_kssVwt_0.1_lfcShrink_genes
do
( gs -q -dNOPAUSE -r300x300 -sDEVICE=tiff24nc -sOutputFile=${i}.tiff ${i}.pdf -c quit

  convert -density 300 ${i}.pdf -quality 90 ${i}.png ) &
done
wait
