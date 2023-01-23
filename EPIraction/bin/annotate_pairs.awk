#!/bin/awk -f

BEGIN { OFS="\t";FS="\t" }

function abs(v) {return v < 0 ? -v : v}

{
    print chrom FS $1 FS ($1+4999) FS $2 FS abs($2-$1) FS $3
    if($1 != $2)
    {
	print chrom FS $2 FS ($2+4999) FS $1 FS abs($1-$2) FS $3
    }
}
