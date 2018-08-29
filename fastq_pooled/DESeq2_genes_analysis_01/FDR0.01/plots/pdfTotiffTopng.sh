#!/bin/bash

for i in MAplot_res_kssVwt_0.01_lfcShrink_genes
do
( convert -density 300 ${i}.pdf -quality 90 ${i}.png ) &
done
wait
