#!/bin/bash

files="input_cms virtequiv pwg-btlgrid.top pwgxgrid.dat pwg-rmngrid.top FlavRegList bornequiv realequivregions-btl pwhg_checklimits"
for fn in $files
do
    rm -rf ./$fn
done
