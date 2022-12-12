import argparse
import pysam
import os
# Adopted by Yang Pan 2020.12.20 (panyang@ucla.edu)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.10.02
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script that extracts splice junctions from a STAR-aligned BAM
# file and annotates each splice junction with the number of uniquely
# mapped reads that support that splice junction. Supporting reads for 
# annotated splice junctions are by default only required to have a 
# minimum overhang of 1 bp. Supporting reads for unannotated and
# canonical splice junctions (GT-AG/CT-AC; GC-AG/CT-GC; AT-AC/GT-AT) are
# by default only required to have a minimum overhang of 8 bp. Supporting
# reads for unannotated and non-canonical splice junctions however are 
# by default required to have a minimum overhang of 10 bp. The resulting
# output file is a TSV with two fields: (1) a splice junction ID with 
# format [chr#]:[start]:[end], (2) number of uniquely mapped reads 
# supporting the splice junction. All genomic coordinates used in 
# describing each splice junction are 1-based. Make sure that the BAM 
# file and genome FASTA file have been indexed. 

# UPDATE (2020.10.02): Users have the option of strictly filtering
# for reads of a specified size using the -r <read_length> option

# Usage:
#   python extract_SJ.py -i /path/to/BAM/file \
#       -g /path/to/annotation/GTF/file \
#       -a [minimum overhang length for annotated SJs, default: 1]
#       -c [minimum overhang length for unannotated canonical SJs, default: 8] \
#       -n [minimum overhang length for unannotated non-canonical SJs, default: 10] \
#       -r [length of reads to keep when counting junction reads] \
#       -f /path/to/genome/fasta/file \
#       -o /path/to/output/file

# Dependencies:
#   * argparse
#   * pysam

def get_introns(blocks):
    introns = []
    
    # N blocks are associated with N-1 introns
    for i in range(len(blocks)-1):
        intronStart = blocks[i][1] + 1
        intronEnd = blocks[i+1][0]
        introns += [(intronStart, intronEnd)]
    
    return introns

def isCanonical(dn1, dn2):
    # Construct tuple representing SJ dinucleotides
    sjDN = (dn1, dn2)
    
    # Establish canonical SJs
    canonicalSJ = [('GT', 'AG'), ('CT', 'AC'), ('GC', 'AG'),
        ('CT', 'GC'), ('AT', 'AC'), ('GT', 'AT')]
    
    return sjDN in canonicalSJ

def get_threshold(chr, introns, genome, annoSJ, minOverhang, minOverhangC, minOverhangNC):
    isAnno = []
    isCanon = []
    
    # Iterate through each intron
    for i in range(len(introns)):
        # Construct splice junction ID
        sjInfo = [chr, introns[i][0], introns[i][1]]
        sjID = ':'.join(map(str, sjInfo))
        
        # Retrieve dinucleotides for splice junction
        dn1 = genome.fetch(chr, introns[i][0]-1, introns[i][0]+1)
        dn2 = genome.fetch(chr, introns[i][1]-2, introns[i][1])
        
        isAnno.append(sjID in annoSJ)
        isCanon.append(isCanonical(dn1, dn2))
        
    if all(isAnno):
        threshold = minOverhang
    elif all(isCanon):
        threshold = minOverhangC
    else:
        threshold = minOverhangNC

    return threshold
    
def update_SJdb(read, sjDB, annoSJ, genome, minOverhang, minOverhangC, minOverhangNC):
    # Get chromosome for read
    chr = read.reference_name
    
    # Extract coordinates of blocks for each read
    blocks = read.get_blocks()
    introns = get_introns(blocks)

    # Determine threshold for read anchor lengths
    threshold = get_threshold(chr, introns, genome, annoSJ, minOverhang, minOverhangC, minOverhangNC)
    
    # Compute anchor lengths
    leftAnchorLen = blocks[0][1] - blocks[0][0]
    rightAnchorLen = blocks[len(introns)][1] - blocks[len(introns)][0]
    
    # Only keep reads that satisfy appropriate anchor length
    # threshold
    if min(leftAnchorLen, rightAnchorLen) >= threshold:
        # Iterate through introns and check if SJ is canonical
        for i in range(len(introns)):
            # Construct splice junction ID
            sjInfo = [chr, introns[i][0], introns[i][1]]
            sjID = ':'.join(map(str, sjInfo))
            
            # Update sjDB with SJ            
            # Check if SJ exists in dictionary
            if sjID in sjDB:
                # Increment current read count by 1
                sjDB[sjID] += 1
            else:
                # Add sjID to sjDB with one read
                sjDB[sjID] = 1

    # Return the updated dictionary
    return sjDB

def get_transcript_ID(infoString):
    infoStringArr = infoString.split(';')
    # Get position in the INFO string that has the transcript ID
    idx = [i for i, s in enumerate(infoStringArr) if 'transcript_id' in s][0]
    return infoStringArr[idx].split('"')[1]

def build_anno_SJdb(gtfPath, chrList):
    annoSJ = {}
    
    # Create dictionary object to store exons of a transcript 
    exons = {}
    # Create dictionary object to store chr of a transcript
    chrDict = {}
    
    # Read through every line of the gtfPath
    with open(gtfPath) as f:
        for line in f:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                if info[0] in chrList and info[2] == 'exon':
                    exonStart = int(info[3])
                    exonEnd = int(info[4])
                    transcriptID = get_transcript_ID(info[8])
                    # Check if transcript is included in exons
                    if transcriptID in exons:
                        exons[transcriptID] += [(exonStart, exonEnd)]
                    else:
                        exons[transcriptID] = [(exonStart, exonEnd)]
                    
                    # Check if transcript is included in chrDict
                    if transcriptID not in chrDict:
                        chrDict[transcriptID] = info[0]

    # Extract splice junctions from exons
    for transcriptID in exons:
        chr = chrDict[transcriptID]
        
        # Sort the tuple of exons
        exonList = sorted(exons[transcriptID])
        
        # Build annoSJ
        for i in range(len(exonList)-1):
            sjID = ':'.join(map(str,[chr, exonList[i][1]+1, exonList[i+1][0]-1]))
            if sjID not in annoSJ:
                annoSJ[sjID] = 0

    return annoSJ

def write_output(sjDB, outfile):
    f = open(outfile, 'w')
    
    # Iterate through all sjIDs in sjDB
    for sjID in list(sjDB.keys()):
        f.write('\t'.join(map(str,[sjID, sjDB[sjID]]))+'\n')

    f.close()

def check_cigar(cigarString):
    # Only keep read if CIGAR string only contains 'N' and 'M'
    return set([i for i in cigarString if i.isalpha()]) == set(['N', 'M'])

def main(args):

    # Parse command-line arguments
    bamPath, gtfPath, fastaPath, outfile = args.bam_path, args.gtf, args.genome_fasta, args.outdir
    minOverhang, minOverhangC, minOverhangNC, filterRL = int(args.minimum_overhang_length_annotated), int(args.minimum_overhang_length_unannotated_canonical), int(args.minimum_overhang_length_unannotated_noncanonical), int(args.read_length)
    
    if os.path.exists(bamPath+'.bai') == False:
        pysam.index(bamPath)
    # Create list of target chromosomes
    chrList = ['chrX', 'chrY']
    for i in range(22):
        chrList.append('chr' + str(i+1))
    
    # Initialize splice junction database (dictionary object)
    sjDB = {}
    
    # Create dictionary object to store all annotated SJs in the GTF
    annoSJ = build_anno_SJdb(gtfPath, chrList)

    # Open the BAM file as a SAM
    samfile = pysam.AlignmentFile(bamPath, 'rb')
    
    # Open FASTA file
    genome = pysam.Fastafile(fastaPath)
    
    # Iterate through the reads in the SAM file and update the 
    # sjDB object
    for read in samfile.fetch():
        # Check if read satisfies the following:
        #   * is proper pair
        #   * is uniquely mapped by STAR
        #   * maps to chromosomes 1-22, X, Y
        #   * CIGAR string only contains M and N
        if (read.is_proper_pair
            and read.mapping_quality == 255
            and read.reference_name in chrList
            and check_cigar(read.cigarstring)):
                # Check if read has the specified read length
                # if provided
                if filterRL == -1:
                    # Update sjDB with given read
                    sjDB = update_SJdb(read, sjDB, annoSJ, genome, 
                        minOverhang, minOverhangC, minOverhangNC)
                else:
                    if read.query_length == filterRL:
                        # Query read has expected length
                        sjDB = update_SJdb(read, sjDB, annoSJ, genome, 
                            minOverhang, minOverhangC, minOverhangNC)
                   
    # Print out sjDB to outfile
    write_output(sjDB, outfile)

if __name__ == '__main__':
    main()
