import pandas as pd
import pyranges as pr
import pybedtools as py
import numpy as np
import argparse
import sys
from pyfaidx import Fasta

def parse_inputs(args_list):
    parser = argparse.ArgumentParser(prog='PROG', usage='Find Danio rerio safe harbors')

    parser.add_argument('-chro', required=True, help='chromosome lengths')
    parser.add_argument('-genes', help='genes annotation file GTF')
    parser.add_argument('-onco', help='List of oncogene names TXT tab')
    parser.add_argument('-enh', help='enhancer coordinates BED')
    parser.add_argument('-cent', help='centromere coordinates GTF')
    parser.add_argument('-gap', help='Gap coordinates BED')
    parser.add_argument('-lnc', help='lncRNA coordinates BED')
    parser.add_argument('-mi', help='miRNA coordinates BED')
    parser.add_argument('-t', help='tRNA coordinates BED')
    parser.add_argument('-rm', help='RepeatMasker coordinates')
    parser.add_argument('-f', help='primary assembly file')


    args = parser.parse_args(args_list)

    return args

#-----------------DISTANCE TABLE--------------------
dist_data = {'Type': ["gene", "oncogene", "microrna", "trna", "lncrna", "enhancer", "centromere", "gap"], 'Distance': [5000, 300000, 300000, 150000, 150000, 20000, 300000, 300000]}
dist_df = pd.DataFrame(data=dist_data).set_index('Type')

#-----------------CHROMOSOME LENGTHS--------------------
def process_chromosome_lengths(filepath):
    print(f"{' PROCESSING CHROMOSOME LENGTHS ':=^60}")
    chrom_len = pd.read_csv(filepath, sep='\t', header=None, names=['Chromosome', 'Lengths'])

    chrom_len['Chromosome'] = chrom_len['Chromosome'].str.strip('chr')
    
    #filter out scaffolds
    chrom_len['Chromosome'] = pd.to_numeric(chrom_len['Chromosome'], errors='coerce') #converts to int if string is number and NaN if not
    chrom_len = chrom_len.dropna(subset=['Chromosome']) #drops all NaN

    #sort chromosomes
    chrom_len = chrom_len.sort_values(by='Chromosome').reset_index(drop=True)
    chrom_len['Chromosome'] = chrom_len['Chromosome'].astype(int).astype(str)

    #set chromosomes as index
    chrom_len = chrom_len.set_index('Chromosome')
    print(f"{' DONE ':=^60}")
    return chrom_len

#-----------------GENES--------------------
def process_genes(filepath):
    print(f"{' PROCESSING GENES ':=^60}")
    #read in gtf file
    gtf = pr.read_gtf(filepath)

    #set as df
    genes_df = gtf.df
    genes_df = genes_df[['Chromosome', 'Feature', 'Start', 'End', 'gene_name']]
    
    #filter out scaffolds
    genes_df['Chromosome'] = pd.to_numeric(genes_df['Chromosome'], errors='coerce') #converts to int if string is number and NaN if not
    genes_df = genes_df.dropna(subset=['Chromosome']) #drops all NaN

    #filter by genes
    genes_df = genes_df[genes_df['Feature'] == 'gene'].reset_index(drop=True)

    #convert to bed format 
    genes_df['Start'] = genes_df['Start'] - 1

    #make the chromosome, start and end columns values numeric
    genes_df[['Chromosome', 'Start', 'End']] = genes_df[['Chromosome', 'Start', 'End']].astype(int)
    genes_df = genes_df.sort_values(by=['Chromosome', 'Start']).reset_index(drop=True)


    print(f"{' DONE ':=^60}")

    return genes_df

#-----------------ONCOGENES--------------------
def process_oncogenes(onco_list, genes_df):
    print(f"{' PROCESSING ONCOGENES ':=^60}")
    
    #read in onco gene names list
    onco_names_df = pd.read_csv(onco_list, sep='\t', header=None, names=['gene_name'])

    #merge with genes on gene_name
    merged = pd.merge(onco_names_df, genes_df, how='inner', on='gene_name')

    #remove gene_name and add feature column
    merged = merged.drop(columns='gene_name')
    merged['Feature'] = 'onco' + merged['Feature']
    
    #sort values
    merged = merged.sort_values(by=['Chromosome', 'Start']).reset_index(drop=True)
    
    print(f"{' DONE ':=^60}")
    
    return merged

#-----------------TRNA, LNCRNA, MIRNA, CENTROMERES, GAPS, ENHANCERS--------------------
def process_others(filepath, feature, columns, d, format):
    print(f"{f' PROCESSING {feature} ':=^60}")
    
    #read in df
    df = pd.read_csv(filepath, sep=d, header=None, usecols=columns, names=['Chromosome', 'Start', 'End'])
    
    #create feature column
    df.insert(1, column='Feature', value=feature)

    #assign column names
    df['Chromosome'] = df['Chromosome'].astype(str)
    
    #filter out the alternates and unknowns
    df = df[~df['Chromosome'].str.contains('_alt|chrUn|_fix')]

    #strip "chr"
    df['Chromosome'] = df['Chromosome'].str.strip('chr')

    #filter out scaffolds
    df['Chromosome'] = pd.to_numeric(df['Chromosome'], errors='coerce') #converts to int if string is number and NaN if not
    df = df.dropna(subset=['Chromosome']) #drops all NaN

    #change datatype of start and end columns
    df['Start'] = df['Start'].astype(int)
    df['End'] = df['End'].astype(int)
    df['Chromosome'] = df['Chromosome'].astype(int)

    #convert to BED format if GTF
    if format == 'GTF':
        df['Start'] = df['Start'] - 1

    #sort by chromo and start col
    df = df.sort_values(by=['Chromosome', 'Start']).reset_index(drop=True)

    print(f"{' DONE ':=^60}")

    return df

#-----------------ADD FLANKS, ADJUST ENDS, MERGE, FIND GSH --------------------
def find_gsh(list_dfs, chromo_lens_df, dist_df, chromo_file):
    temp_df = []
    
    for df in list_dfs:
        df = df.copy()
        
        #check if feature is in distance table
        feature = df.iloc[0, 1]
        if feature in dist_df.index:
            flank = dist_df.loc[feature, 'Distance']

            print(f"{f' ADDING FLANKS AND ADJUSTING ENDS {feature} ':=^60}")
            print('Feature:', feature, 'Rows:', df.shape[0], 'Columns:', df.shape[1])
            
            #add flanks
            df['Start'] = df['Start'] - flank
            df['End'] = df['End'] + flank

            #change start to 0 is negative 
            df.loc[df['Start'] < 0, 'Start'] = 0

            #adjust ends if exceeds chromosome lengths
            chromo_lens_df.index = chromo_lens_df.index.astype(int)
            df = df.merge(chromo_lens_df, left_on='Chromosome', right_index=True)
            df['End'] = df[['End', 'Lengths']].min(axis=1)
            df = df.drop(columns='Lengths')
            chromo_lens_df.index = chromo_lens_df.index.astype(str)

            temp_df.append(df)

            print(f"{f' DONE WITH {feature} ':=^60}")
        else:
            print("ERROR Feature: ", feature, 'not found!')

    #concat add dfs together
    regions_to_avoid = pd.concat(temp_df, ignore_index=True)

    #remove columns that are no longer needed
    regions_to_avoid = regions_to_avoid.drop(columns=['Feature', 'gene_name'])
    regions_to_avoid['Chromosome'] = regions_to_avoid['Chromosome'].astype(str)

    #add 'chr' to chromosome number for pybedtools
    regions_to_avoid['Chromosome'] = 'chr' + regions_to_avoid['Chromosome']

    print(f"{' MERGING ALL REGIONS TO AVOID AND FINDING SAFEHARBORS':=^60}")
    
    #create pybedtool 
    concat_bed = py.BedTool.from_dataframe(regions_to_avoid)

    #sort based on chromosome lengths
    sorted_bed = concat_bed.sort(g=chromo_file)

    #merge overlaps
    merged = sorted_bed.merge()

    #find complement
    safe_harbors = merged.complement(g=chromo_file)
    
    #save to dataframe
    safe_harbors_df = safe_harbors.to_dataframe(names=['Chromosome', 'Start', 'End'])
    
    #add length of safeharbors
    safe_harbors_df['Size in bp'] = safe_harbors_df['End'] - safe_harbors_df['Start']

    #save to tsv file
    safe_harbors_df.to_csv('safeharbors.tsv', sep='\t', index=False)

    print(f"{f' DONE: Found {safe_harbors_df.shape[0]} safeharbors - Saved to safeharbors.tsv ':=^60}")
    
    return safe_harbors, safe_harbors_df

#-----------------EXTRACTING GSH SEQUENCES--------------------
def extract_seqs(safe_harbors, fasta, filename):
    print(f"{' EXTRACTING SAFEHARBOR SEQUENCES ':=^60}")
    
    genome = Fasta(fasta)
    with open(filename, 'w') as f:
        for i, row in safe_harbors.iterrows():
            key = row['Chromosome'].lstrip('chr')
            start = row['Start']
            end = row['End']
            length = end - start

            gsh_seq = genome[key][start:end]
            f.write(f'>chr{gsh_seq.name}:{gsh_seq.start}-{gsh_seq.end}|{length}bp\n')
            f.write(gsh_seq.seq + '\n')

    print(f"{f' DONE: safeharbor sequences saved to {filename} ':=^60}")

#-----------------SUBTRACT RM REGIONS--------------------
def subtract_rm(safe_harbors, rmsk_path, chromo_file):
    print(f"{' SUBTRACTING REPEAT MASKER COORDINATES ':=^60}")
    
    rmsk = py.BedTool(rmsk_path)

    #subtract RM data
    subtracted = safe_harbors.subtract(rmsk).sort(g=chromo_file)
    
    #save as dataframe
    df = subtracted.to_dataframe(names=['Chromosome', 'Start', 'End'])
    
    #filter results to regions larger than 100 bp
    df = df[df['End'] - df['Start'] > 100].reset_index(drop=True)
    df['Size in bp'] = df['End'] - df['Start']
    df.to_csv('safe_harbors_rm.tsv', sep='\t', index=False)

    print(f"{' DONE: safeharbors with subtracted RepeatMasker coordinates saved to safeharbors_rm.tsv ':=^60}")
    
    return subtracted, df

def main():
    #Parses inputs from command-line.
    args = parse_inputs(sys.argv[1:])
    
    #chromosome lengths
    chromo_lens = process_chromosome_lengths(args.chro)

    #genes
    genes = process_genes(args.genes)

    #oncogenes
    oncogenes = process_oncogenes(args.onco, genes)

    #enhancers
    enhancers = process_others(args.enh, "enhancer", [0, 1, 2], '\t', "BED")
    
    #centromeres
    centromeres = process_others(args.cent, 'centromere', [0, 1, 2], '\s+', 'GTF')

    #gaps
    gaps = process_others(args.gap, 'gap', [0, 1, 2], '\s+', 'BED')

    #lncRNAs
    lncrnas = process_others(args.lnc, "lncrna", [0, 1, 2], '\t', 'BED')

    #miRNAs
    mirnas = process_others(args.mi, "microrna", [0, 1, 2], '\t', 'BED')

    #tRNAs
    trnas = process_others(args.t, "trna", [0, 3, 4], '\t', 'GTF')

    #safe harbors
    gsh = find_gsh([genes, enhancers, oncogenes, trnas, mirnas, lncrnas, centromeres, gaps], chromo_lens, dist_df, args.chro)

    #sequences
    seqs = extract_seqs(gsh[1], args.f, 'safeharbors.fasta')

    #rm safe harbors
    gsh_rm = subtract_rm(gsh[0], args.rm, args.chro)

    #rm safe harbor sequences
    gsh_rm_seqs = extract_seqs(gsh_rm[1], args.f, 'safeharbors_rm.fasta')

if __name__ == "__main__":
	main()