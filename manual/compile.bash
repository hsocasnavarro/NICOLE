#!/bin/bash
latex manual
latex manual
dvips manual -o
ps2pdf manual.ps manual.pdf
cp manual.pdf ../
#rm *~
find ../ -name "*~" -exec rm {} \;