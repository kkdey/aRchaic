from pyfaidx import Fasta
import csv
import pysam
import argparse as arg
import random

def find_substitutions_wtags(aligned_pairs):
    """
    For every mutation, the mutation position
    along the read (mpos) and the position of the mutation in the genome (posg)
    """
    mpos = ()
    posg = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                mpos = mpos + (aligned_pairs[i][0],)
                posg = posg + (aligned_pairs[i][1],)
                
    return((mpos, posg))

def find_substitutions(read, reference):
    # chrom is 0-based
    chrom = read.reference_id
    ref_pos = read.get_reference_positions()
    seq = read.query_alignment_sequence

    if (ref_pos is None or None in ref_pos or len(ref_pos) is 0):
        return((), ())

    start = ref_pos[0]
    end = ref_pos[-1]
    refs = reference[chrom][start:(end+1)]

    # check if it's an indel
    if ((end-start+1) is not len(seq) or len(seq) is not len(refs) or len(refs) is not len(ref_pos)):
        return((), ())

    mpos = [i for i in range(len(refs)) if refs[i] != seq[i]]
    posg = [ref_pos[pos] for pos in mpos]

    return((mpos, posg))

def get_bases_strandbreaks(read, reference):
    # chrom is 0_based
    chrom = read.reference_id
    ref_pos = read.get_reference_positions()
    leftbreak = reference[chrom][ref_pos[0]-1]
    rightbreak = reference[chrom][ref_pos[-1]+1]
    return((leftbreak, rightbreak))

if __name__ == '__main__':
    parser = arg.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="bam file, the bam file must be indexed")
    parser.add_argument("-f", "--fasta", required=True, help = "reference file, the reference file must be indexed")
    parser.add_argument("-o", "--out", required=True, help="out file")
    parser.add_argument("--dont-use-tags", help = "Don't use the MD and NM tags? If the tags are computed incorrectly from the alignment, you can use this option (slower option)", default = False, action = 'store_true', dest = "dont_use_tags")
    parser.add_argument("-n", "--nflanking", required=False, default = 1, help="number of flanking base-pairs around the mismatch to record, must be an integer >= 1", dest = "nflanking", type=int)
    parser.add_argument("-mq", "--mapquality", required=False, default = 20, help="min. mapping quality", dest = "mq", type = int)
    parser.add_argument("-bq", "--basequality", required=False, default = 30, help="min. base quality", dest = "bq", type = int)
    parser.add_argument("-ss", "--num-subsample", required=False, default = 20000, help="number of reads to subsample", dest = "n_subsample", type = int)
    parser.add_argument('-seed', '--seed', required = False, default = 100, help = "seed value", dest = 'seed', type = int)
    
    args = parser.parse_args()

    ## ../data/bam/Lindo2016/T004_all_chr.bam
    bam = pysam.AlignmentFile(args.bam, "rb")
    ## "/project/jnovembre/data/external_public/reference_genomes/hs37d5.fa"
    fastafile = Fasta(args.fasta, as_raw = True)

    mismatches = []
    
    ## Here, we borrow some of John Blishak's code to sample reads
    ## https://github.com/jdblischak/singleCellSeq/blob/master/code/subsample-bam.py

    num = min(args.n_subsample, bam.mapped)

    indices = list(range(bam.mapped))
    random.seed(args.seed)
    random.shuffle(indices)
    subbed = indices[0:num]
    subbed.sort()
    reads = list()

    index = 0
    list_index = 0
    for read in bam.fetch(until_eof = True):
        if (read.is_unmapped):
            continue

        if (read.mapping_quality < args.mq or read.is_duplicate):
            continue
        
        if (not args.dont_use_tags) and read.get_tag('NM') == 0:
            continue

        if index == subbed[list_index]:
            reads.append(read)
            list_index += 1
            if list_index == num:
                break
        index += 1

    
    for r in range(len(reads)):

        read = reads[r]

        seq = read.query_sequence

        if (args.dont_use_tags):
            (mutPos, posg) = find_substitutions(read, fastafile)
        else:
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            (mutPos, posg) = find_substitutions_wtags(aligned_pairs)

            
        baseq = read.query_qualities

        if (read.is_reverse):
            strando = '-'
        else:
            strando = '+'

        try:
            (leftbreak, rightbreak) = get_bases_strandbreaks(read, fastafile)
        except IndexError:
            leftbreak = 'N'
            rightbreak = 'N'
                
        for i in range(len(posg)):
            pos = mutPos[i]
            mut = seq[pos]
        
            if (baseq[pos] < args.bq):
                continue

            # we don't count mutation in soft clipped areas
            if (pos < read.qstart or pos > read.qend or mut == 'N'):
                continue

            start = posg[i]-args.nflanking
            end = posg[i]+args.nflanking
            ref = fastafile[read.reference_id][start:(end+1)]
            patt = (ref[0:(args.nflanking+1)]).upper() + '->' + mut.upper() + (ref[(args.nflanking+1):(2*args.nflanking+1)]).upper()

            mutStart = pos - read.query_alignment_start
            mutEnd = read.qlen - (pos + 1)

            val = (patt, mutStart, mutEnd, leftbreak.upper(), rightbreak.upper(), strando, r)
            mismatches.append(val)
          
    # write to file
    with open(args.out, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for val in mismatches:
            writer.writerow([val[0], val[1], val[2], val[3], val[4], val[5], val[6]])


