import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from scipy.cluster import hierarchy
import io
import csv

def pars():

    parser = argparse.ArgumentParser(prog="cluster.py")
    parser.add_argument("-v","--vcf",dest="vcf",help="Filtered VCF File")
    parser.add_argument("-y","--height",dest="height",help="Height to cut dendrogram",type=int,default=None)
    parser.add_argument("-n","--num_groups",dest="groups",help="Number of groups desired",type=int,default=None)
    parser.add_argument("-o","--outfolder",dest="outfolder",help="Folder name for output",default='recomboutput')
    args = parser.parse_args()

    outfolder = args.outfolder
    outfolder = outfolder+"/0.cluster/"

    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    return args.vcf, outfolder, args.height, args.groups

def read_vcf(path): #from https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def editvcf(vcf):

    vcf.drop(vcf.iloc[:, 0:9], inplace = True, axis = 1)    #keep only allele data
    for name in list(vcf.columns):
        if vcf[name].dtypes != 'int64':
            vcf[name] = vcf[name].str[:3]                   #keep only genotype data
    vcf = vcf.replace(['1|1','1/1'],'2')                    #convert allele calls to counts of alt alleles
    vcf = vcf.replace(['0|0','0/0'],'0')
    vcf = vcf.replace(['.|.','./.'],'-1')
    vcf = vcf.replace(['0|1','1|0','1/0','0/1'],'1')
    vcf = vcf.astype(int)                                   #convert to int

    return vcf

def heatmap(vcf,outfolder,vcfname,chroms,pos,minorticks):

    print("\n\tClustering...")

    fig, ax = plt.subplots(figsize=(11,20))
    colours = ['w','b','y','r']
    custom_palette = sns.set_palette(sns.color_palette(colours))
    palette = sns.color_palette(custom_palette, desat=0.65, as_cmap=True)

    before = sns.clustermap(vcf,method='ward',cmap=palette,vmin=-1,vmax=2,cbar_pos=None,xticklabels=False,row_cluster=False,col_cluster=False,dendrogram_ratio=0.0)

    outbefore = os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_before_clustering.png")
    plt.ylabel('Chromosome')
    plt.xlabel('Individual')
    legend_elements = [Patch(facecolor='b', label='Homozygous ref'), Patch(facecolor='r', label='Homozygous alt'), Patch(facecolor='y', label='Heterozygous'), Patch(facecolor='w', edgecolor='gray', label='Missing')]
    plt.legend(handles=legend_elements, ncol=4, bbox_to_anchor=(.5, 1.03),loc='center')
    plt.yticks(pos,chroms)
    plt.tight_layout()
    plt.savefig(outbefore, dpi=600,bbox_inches='tight')
    print("\n\tBefore clustering saved as",outbefore)

    after = sns.clustermap(vcf,method='ward',cmap=palette,vmin=-1,vmax=2,cbar_pos=None,xticklabels=False,row_cluster=False,dendrogram_ratio=0.1)

    outafter = os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_after_clustering.png")
    plt.ylabel('Chromosome')
    plt.xlabel('Individual')
    plt.yticks(pos,chroms)
    plt.tight_layout()
    plt.savefig(outafter, dpi=600,bbox_inches='tight')
    print("\n\tAfter clustering saved as",outafter)

    leaves = after.dendrogram_col.dendrogram['leaves']
    dgram2 = after.dendrogram_col.linkage

    plt.figure(figsize=(25,5))
    dn = hierarchy.dendrogram(dgram2,labels=list(vcf.columns.values),get_leaves=True)
    out_dendrogram = os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_dendrogram.png")
    plt.ylabel('Height')
    plt.xlabel('Individual')
    plt.title('Dendrogram')
    plt.tight_layout()
    plt.savefig(out_dendrogram, dpi=600,bbox_inches='tight')
    print("\n\tDendrogram saved as",out_dendrogram)

    out_linkage = os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_linkage")
    np.save(out_linkage, dgram2)


def sets(vcf,outfolder,vcfname,chroms,pos,ngroups=None,h=None):

    print("\n\tLoading saved data...")
    dgram2 = np.load(os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_linkage.npy"),allow_pickle='TRUE')

    leavesindex = hierarchy.leaves_list(dgram2)
    leaveslist = vcf.columns[leavesindex]
    leavesedit = []
    for name in leaveslist:
        leavesedit.append(name.split('__TCAP',1)[0])

    print('\n\tCutting tree...')

    if h is not None:
        plt.figure(figsize=(10,5))
        dn_cut = hierarchy.dendrogram(dgram2,labels=leavesedit,get_leaves=True,leaf_font_size=2.5)
        plt.axhline(y = h, color = 'r', linestyle = 'dashed')
        out_cut_dendrogram = os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_cut_dendrogram.png")
        plt.ylabel('Height')
        plt.xlabel('Individual')
        plt.tight_layout()
        plt.savefig(out_cut_dendrogram, dpi=600,bbox_inches='tight')
        print("\n\tCut dendrogram saved as",out_cut_dendrogram)

    cut = hierarchy.cut_tree(dgram2,n_clusters=ngroups,height=h)
    cut2 = np.squeeze(cut)

    leaves2 = []
    x = 0
    while x < len(leaveslist):
        leaves2.append([leaveslist[x],cut2[x]])
        x+=1

    leaves2.sort(key=lambda x: x[1])

    #save sets file
    outfile = os.path.join(outfolder,os.path.basename(vcfname)[:-4]+"_sets.txt")
    with open(outfile,'w+') as out_sets:
        writer = csv.writer(out_sets, delimiter="\t")
        for row in leaves2:
            writer.writerow(row)
    out_sets.close()

    print("\n\tSets file saved as",outfile)

def main():

    print("cluster.py:")
    vcfname, outfolder, height, groups = pars()
    print("\n\tInput vcf:",vcfname,"\n\tOutput folder:",outfolder)

    vcf = read_vcf(vcfname)

    illegal = ["Pan_BP_CH","chrUn"]
    for chr in illegal:
        vcf.drop(vcf.loc[vcf['CHROM']==chr].index, inplace=True)
    chroms = vcf['CHROM'].unique()

    lenlist = [0]
    for chr in chroms:
        lenlist.append(len(vcf[vcf['CHROM']==chr]))

    lenlist1 = []
    for x in range(len(lenlist)):
        lenlist1.append(sum(lenlist[:x+1]))

    lenlist2 = []
    x = 1
    while x < len(lenlist1):
        lenlist2.append((lenlist1[x-1]+lenlist1[x])/2)
        x+=1

    if height or groups is not None:
        sets(editvcf(vcf),outfolder,vcfname,chroms,lenlist2,groups,height)
    else:
        heatmap(editvcf(vcf),outfolder,vcfname,chroms,lenlist2,lenlist1)
    print("\ncluster.py finished.\n")

if __name__ == '__main__':
    main()
