from pyfaidx import Fasta
import csv
import pysam
import argparse as arg

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

def find_substitutions(read, chr, reference):
    ref_pos = read.get_reference_positions()
    seq = read.query_alignment_sequence

    if (ref_pos is None or None in ref_pos or len(ref_pos) is 0):
        return((), ())

    start = ref_pos[0]
    end = ref_pos[-1]
    refs = reference[(chr-1)][start:(end+1)]

    # check if it's an indel
    if ((end-start+1) is not len(seq) or len(seq) is not len(refs) or len(refs) is not len(ref_pos)):
        return((), ())

    mpos = [i for i in range(len(refs)) if refs[i] != seq[i]]
    posg = [ref_pos[pos] for pos in mpos]

    return((mpos, posg))

def get_bases_strandbreaks(read, chr, reference):
    ref_pos = read.get_reference_positions()
    leftbreak = reference[(chr-1)][ref_pos[0]-1]
    rightbreak = reference[(chr-1)][ref_pos[-1]+1]
    return((leftbreak, rightbreak))

if __name__ == '__main__':
    parser = arg.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="bam file, the bam file must be indexed")
    parser.add_argument("-f", "--fasta", required=True, help = "reference file, the reference file must be indexed")
    parser.add_argument("-o", "--out", required=True, help="out file")
    parser.add_argument("--add-chr", help = "Does your reference use the chr prefix? For example, different versions of the references either use chr1 or 1 to designate chromosome 1, you can find out by running samtools idxstats <your bamfile> | head -1 ", default = False, action = 'store_true', dest = "add_chr")
    parser.add_argument("--dont-use-tags", help = "don't use the MD and NM tags? If the tags are computed incorrectly from the alignment, you can use this option (slower option)", default = False, action = 'store_true', dest = "dont_use_tags")
    parser.add_argument("-n", "--nflanking", required=False, default = 1, help="number of flanking base-pairs around the mismatch to record, must be an integer >= 1", dest = "nflanking", type=int)
    parser.add_argument("-mq", "--mapquality", required=False, default = 20, help="min. mapping quality", dest = "mq", type = int)
    parser.add_argument("-bq", "--basequality", required=False, default = 30, help="min. base quality", dest = "bq", type = int)
    #parser.add_argument("-id", "--add-read-id", required=False, default = False, help="add a unique read id?", action = 'store_true', dest = "add_id")
    
    args = parser.parse_args()

    ## ../data/bam/Lindo2016/T004_all_chr.bam
    samfile = pysam.AlignmentFile(args.bam, "rb")
    ## "/project/jnovembre/data/external_public/reference_genomes/hs37d5.fa"
    fastafile = Fasta(args.fasta, as_raw = True)

    mismatches = []
    
    chrs = [i for i in range(1,23)]

    read_cnt = 0
    for chr in chrs:
        
        for read in samfile.fetch(('chr' + str(chr)) if args.add_chr else str(chr)):
            if (read.mapping_quality < args.mq or read.is_duplicate):
                continue

            if (not args.dont_use_tags) and read.get_tag('NM') == 0:
                continue

            read_cnt += 1

            seq = read.query_sequence

            if (args.dont_use_tags):
                (mutPos, posg) = find_substitutions(read, chr, fastafile)
            else:
                aligned_pairs = read.get_aligned_pairs(with_seq=True)
                (mutPos, posg) = find_substitutions_wtags(aligned_pairs)

            
            baseq = read.query_qualities

            if (read.is_reverse):
                strando = '-'
            else:
                strando = '+'

            try:
                (leftbreak, rightbreak) = get_bases_strandbreaks(read, chr, fastafile)
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
                ref = fastafile[(chr-1)][start:(end+1)]
                patt = (ref[0:(args.nflanking+1)]).upper() + '->' + mut.upper() + (ref[(args.nflanking+1):(2*args.nflanking+1)]).upper()

                mutStart = pos - read.query_alignment_start
                mutEnd = read.qlen - (pos + 1)

                read_name = read_cnt
                
                val = (patt, mutStart, mutEnd, leftbreak.upper(), rightbreak.upper(), strando, read_name)
                mismatches.append(val)
          
    # write to file
    with open(args.out, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for val in mismatches:
            writer.writerow([val[0], val[1], val[2], val[3], val[4], val[5], val[6]])


