#!/bin/bash


#flist="/home/rj18633/scratch/gp.map/data/besd.formatting/test.out/v1/v3/flist/test.flist"
#besd_out="/home/rj18633/scratch/gp.map/data/besd.formatting/test.out/v1/v3/besd/test_pipeline"

flist="$FLIST"
besd_out="$BESD_OUT"

smr --eqtl-flist "$flist" --make-besd --out "$besd_out"
