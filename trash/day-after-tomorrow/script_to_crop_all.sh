#!/bin/bash

for pdf_file in *pdf
do
	pdfcrop $pdf_file $pdf_file
done

cp *pdf /home/lluis/paper_overleaf/low_order_fullBiot/Fig/

