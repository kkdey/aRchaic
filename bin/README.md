How to make the csv file to input in R

(1) index the file with samtools

`samtools index example.bam`

(2) Run the python scipt with default options

`python generate_summary_bams.py -b example.bam -f /path/to/reference/hs37d5\
.fa -o example.csv`


To view the possible options, please type

`python --help generate_summary_bams.py`

