#!/bin/bash


flist="$FLIST"
besd_out="$BESD_OUT"

smr --eqtl-flist "$flist" --make-besd --out "$besd_out"
