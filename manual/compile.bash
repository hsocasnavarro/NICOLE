#!/bin/bash
pdflatex manual
cp manual.pdf ../
#rm *~
find ../ -name "*~" -exec rm {} \;
