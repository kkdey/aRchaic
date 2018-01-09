import pysam
import sys

bamfile = sys.argv[1]
sam = pysam.AlignmentFile(bamfile, "rb")

for read in sam.fetch():
    if read.is_paired:
        print(read.isize, read.qlen)

        
    
