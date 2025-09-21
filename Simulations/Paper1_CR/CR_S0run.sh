#!/bin/sh
# On scripts folder
module load r #/4.1.3
Rscript $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} #> Slurm/output_$2$3$4$5$6$7$8$9$10.txt
