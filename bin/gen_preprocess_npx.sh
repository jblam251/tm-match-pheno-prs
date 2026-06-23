#!/usr/bin/bash

## script will import plink format npx and generate various 
## pre-processed versions of the npx for downstream testing

## INPUT 
##  - NPX in plink format

## OUTPUTS
##  - Impute    + PC0  + No.RINT
##  - Impute    + PC0  + RINT
##  - Impute    + PCn  + No.RINT
##  - Impute    + PCn  + RINT
##  - No.impute + PC0  + No.RINT
##  - No.impute + PC0  + RINT
