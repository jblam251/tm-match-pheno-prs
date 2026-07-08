#!/bin/bash

## USAGE
## $ subset.prs.sh [PRS] [IDMAP] [OUT]

PRS=$1
IDMAP=$2
OUT=$3

zcat $PRS | head -n1 > $OUT
zcat $PRS | grep -f <(cut -f1 $IDMAP) >> $OUT
