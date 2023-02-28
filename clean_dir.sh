#!/bin/bash

rm -r opt
rm -r results
shopt -s extglob
cd resources && rm -r !(adaptive_gls.csv)
