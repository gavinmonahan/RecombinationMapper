import allel
import numpy as np
import sys
import argparse
import csv
import os
import pandas as pd

def pars():
    parser = argparse.ArgumentParser(prog="allele_blocks.py")
    parser.add_argument("-v","--vcf",dest="vcf",help="Filtered VCF File")
    parser.add_argument("-o","--outfolder",dest="outfolder",help="Folder name for output",default='recomboutput')
    parser.add_argument("-p","--phased",dest="phased",help="If genotypes are ALL phased", default=False, action='store_true')
    args = parser.parse_args()

    outfolder = args.outfolder
    outfolder = outfolder+"/1.allele_blocks/"

    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    return args.vcf,outfolder,args.phased

def getcall(call):
    if call[0] == call[1]:
        if call[0] == 0:
            return 0
        if call[0] == 1:
            return 1
        if call[0] == -1:
            return -1
    elif -1 in call:
        if call[0] == -1:
            return call[1]
        else:
            return call[0]
    elif call[0] != call[1]:
        return 2

def blocks(chr, lis, window, sample, positions):
    blocklist = []
    x = 0
    y = 1
    if lis[-1]==0:
        lis.append(1)
    else:
        lis.append(0)
    while x+y < len(lis):
        if lis[x] == lis[x+y]:
            blockstart = x
            y+=1
            continue
        elif lis[x+y] == 2 or lis[x+y] == -1:
            blockstart = x
            y+=1
            continue
        else:
            blockend = x+y-1
            if y >= window:
                blocklist.append([chr,sample,lis[x],blockstart,blockend,positions[blockstart],positions[blockend],positions[blockend]-positions[blockstart],y,window])
            x = x+y
            y = 1
            continue
    return blocklist

def main():
    print("allele_blocks.py:")
    invcf,outfolder,phased = pars()
    print("\n\tInput VCF: ",invcf)
    print("\tOutput folder: ",outfolder)

    header = allel.read_vcf_headers(invcf)
    GT = 0
    for h in header:
        for h2 in h:
            if "Genotype" in h2:
                GT +=1
    if GT == 0:
        print("\tHeader error: check if Genotype format tag is in header.")
        sys.exit(1)

    if phased == True:
        print('\tAssuming phased genotypes.')
    else:
        print('\tAssuming unphased genotypes.')

    print("\n\tFinding blocks:")

    inputvcf = allel.read_vcf(invcf, fields=['samples','variants/CHROM','variants/POS'])
    allchrs = inputvcf['variants/CHROM']
    chroms = np.unique(allchrs)
    illegal=["Pan_BP_CH","chrUn"]
    chrs=[]
    for chr in chroms:
        if chr in illegal:
            continue
        else:
            chrs.append(chr)

    positions = inputvcf['variants/POS']
    minpos = 0
    maxpos = np.amax(positions)
    samples = inputvcf['samples']

    missing = []

    windowsize = 2
    head = ["Chr","Sample","Allele","Block start","Block end","Start pos","End pos","Length","Num SNPs","Win size"]
    outfile = os.path.join(outfolder,os.path.basename(invcf)[:-4]+"_blocks.txt")

    with open(outfile,'w+') as out_co:
        writer = csv.writer(out_co, delimiter="\t")
        writer.writerow(head)

        for chr in chrs:
            readin = allel.read_vcf(invcf, region='{}:{}-{}'.format(chr,minpos,maxpos))
            variants = allel.GenotypeArray(readin['calldata/GT'])
            pos = (readin['variants/POS'])
            miss = variants.count_missing(axis=1)

            for missed in miss:
                for pos3 in pos:
                    missing.append([chr,pos3,missed])

            for col in range(len(variants[0,:])):
                calllist = []
                calllist2 = []
                sample = samples[col]

                for call in variants[:,col]:
                    if phased == True:
                        calllist.append(call[0])
                        calllist2.append(call[1])
                    elif phased == False:
                        calllist.append(getcall(call))

                for row in blocks(chr,calllist,2,sample,pos):
                    writer.writerow(row)
                if phased == True:
                    for row in blocks(chr,calllist2,2,sample+"_2",pos):
                        writer.writerow(row)

            print("\t",chr, " done")
    print("\n\tBlocks saved as", outfile)

    outfile2 = os.path.join(outfolder,os.path.basename(invcf)[:-4]+"_missing.txt")
    missingdf = pd.DataFrame(missing,columns = ['Chr','Pos','Missingness'])
    missingdf.to_csv(outfile2,sep='\t',index=False)
    print("\n\tMissingness saved as", outfile2)

    print("\nallele_blocks.py finished.\n")

if __name__ == '__main__':
    main()
