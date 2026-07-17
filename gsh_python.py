import pandas as pd
import pyranges as pr
import pybedtools as py
import numpy as np
import argparse
import sys
from pyfaidx import Fasta
import subprocess
from pathlib import Path

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
    parser.add_argument('-len', help='minmum length of genomic safeharbor')


    args = parser.parse_args(args_list)

    return args

#-----------------DISTANCE TABLE--------------------
dist_data = {'Type': ["gene", "oncogene", "microrna", "trna", "lncrna", "enhancer", "centromere", "gap", "rm", "non-sunique"], 
            'Distance': [5000, 300000, 300000, 150000, 150000, 20000, 300000, 300000, 0, 0]}
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

    #save regions to avoid for future use
    regions_to_avoid_df = merged.to_dataframe(names=['Chromosome', 'Start', 'End'])

    #find complement
    safe_harbors = merged.complement(g=chromo_file)

    #save as a dataframe
    safe_harbors_df = safe_harbors.to_dataframe(names=['Chromosome', 'Start', 'End'])
    
    #add LP number column
    safe_harbors_df['LP'] = range(1, len(safe_harbors_df) + 1)

    #reorder columns
    safe_harbors_df = safe_harbors_df[['Chromosome', 'Start', 'End', 'LP']]

    #create output folder
    output_dir = Path('output')
    output_dir.mkdir(parents=True, exist_ok=True)

    #save to tsv file and add length of safeharbors
    safe_harbors_df['Length'] = safe_harbors_df['End'] - safe_harbors_df['Start']
    safe_harbors_df = safe_harbors_df[['Chromosome', 'Start', 'End', 'Length', 'LP']]
    safe_harbors_df.to_csv('output/safeharbors.tsv', sep='\t', index=False)

    print(f"{f' DONE: Found {safe_harbors_df.shape[0]} safeharbor regions - Saved to safeharbors.tsv ':=^60}")
    
    return safe_harbors_df, regions_to_avoid_df

#-----------------EXTRACTING GSH SEQUENCES--------------------
def extract_seqs(safe_harbors, fasta, filename):
    print(f"{' EXTRACTING SAFEHARBOR SEQUENCES ':=^60}")
    
    genome = Fasta(fasta)
    with open(filename, 'w') as f:
        for i, row in safe_harbors.iterrows():
            lp_num = row['LP']
            key = row['Chromosome'].lstrip('chr')
            start = row['Start']
            end = row['End']
            length = row['Length']

            gsh_seq = genome[key][start:end]
            f.write(f'>chr{gsh_seq.name}_{gsh_seq.start}_{gsh_seq.end}_{length}bp_LP{lp_num}\n')
            f.write(gsh_seq.seq + '\n')

    print(f"{f' DONE: safeharbor sequences saved to {filename} ':=^60}")

#-------PROCESS REPEAT MASKER DATA---------
def process_rm(rm_file):
    print(f"{' PROCESSING REPEAT MASKER DATA ':=^60}")

    rm = pd.read_csv(rm_file, sep='\t', usecols=[0, 1, 2], header=None, names=['Chromosome', 'Start', 'End'])
    #filter out the alternates and unknowns
    rm = rm[~rm['Chromosome'].str.contains('_alt|chrUn|_fix')]

    #strip "chr"
    rm['Chromosome'] = rm['Chromosome'].str.strip('chr')

    #filter out scaffolds
    rm['Chromosome'] = pd.to_numeric(rm['Chromosome'], errors='coerce') #converts to int if string is number and NaN if not
    rm = rm.dropna(subset=['Chromosome']) #drops all NaN

    rm['Chromosome'] = rm['Chromosome'].astype(int) 
    rm['Chromosome'] = 'chr' + rm['Chromosome'].astype(str)
    
    print(f"{' DONE ':=^60}")

    return rm

#-----ADD REGIONS TO AVOID---------
def add_regions_to_avoid(regions_to_avoid, more_regions_list, chromo_file, min_len, filename):

    print(f"{' ADDING ADDITIONAL REGIONS TO AVOID ':=^60}")

    all_dfs = [regions_to_avoid] + [df.copy() for df in more_regions_list]

    all_regions = pd.concat(all_dfs, ignore_index=False)

    all_regions_pybed = py.BedTool.from_dataframe(all_regions)

    sorted = all_regions_pybed.sort(g=chromo_file)
    
    merged = sorted.merge()
    
    gsh = merged.complement(g=chromo_file)

    gsh_df = gsh.to_dataframe(names=['Chromosome', 'Start', 'End'])
    gsh_df['Start'] = gsh_df['Start'] + 1
    gsh_df['End'] = gsh_df['End'] - 1

    #filter results to regions larger than 100 bp
    gsh_df = gsh_df[gsh_df['End'] - gsh_df['Start'] > min_len].reset_index(drop=True)
    gsh_df['Length'] = gsh_df['End'] - gsh_df['Start']

    #add LP number column
    gsh_df['LP'] = range(1, len(gsh_df) + 1)

    gsh_df.to_csv(filename, sep='\t', index=False)

    print(f"{f' DONE: Found {gsh_df.shape[0]} safeharbor regions - Saved to {filename} ':=^60}")

    return gsh_df

#---------RUN LOCAL BLAST-------------
def run_blast(gsh_file, e_value, outpath):
    print(f"{f' RUNNING LOCAL BLAST ON {gsh_file} ':=^60}")
    #run blast locally 
    db_path = 'danRer11_db/danRer11' #hardcoded path for BLAST database
    subprocess.run([
        "blastn",
        "-query", gsh_file,
        "-db", db_path,
        "-evalue", str(e_value),
        "-out", outpath,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"], check=True)
    
    print(f"{' DONE ':=^60}")

#-----------FIND NON-UNIQUE REGIONS------------
def find_non_unique_seqs(blast_output):
    print(f"{' FINDING NON-UNIQUE SEQUENCES ON SAFEHARBORS ':=^60}")
    df = pd.read_csv(blast_output, sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle'])
    
    #filter out scaffolds
    df['sseqid'] = pd.to_numeric(df['sseqid'], errors='coerce') #converts to int if string is number and NaN if not
    df = df.dropna(subset=['sseqid']) #drops all NaN
    
    #cast as int datatype
    df = df.astype({'qstart': int, 'qend': int})

    #filter out self hits
    for i, row in df.iterrows():
        qseqid = row['qseqid']
        sub_qseqid = qseqid.split('_')

        chrom = int((sub_qseqid[0]).strip('chr'))
        chrom_s = int(row['sseqid'])

        #if chromosome are equal, check if coordinates are contained within 
        if chrom == chrom_s:
            gstart = int(sub_qseqid[1])
            gend = int(sub_qseqid[2])
            sstart = min(row['sstart'], row['send'])
            send = max(row['sstart'], row['send'])

            if not (gend < sstart or gstart > send):
                df.drop(i, inplace=True)

    #preprocess
    non_unique_regions = df[['qseqid', 'qstart', 'qend']].copy()
    qseqid = non_unique_regions['qseqid'].str.split('_')
    non_unique_regions['LP'] = qseqid.str[4] #landing pad number
    non_unique_regions['qseqid'] = qseqid.str[0].str.strip('chr') #chromosome number
    non_unique_regions = non_unique_regions.astype({'qstart': int, 'qend': int}) #cast as ints

    #shift to global coordinates
    non_unique_regions['gstart'] = qseqid.str[1].astype(int)
    #start coordinates
    non_unique_regions['qstart'] = np.where(
                                            non_unique_regions['qstart'] == 1, #condtion
                                            non_unique_regions ['gstart'], #if true
                                            non_unique_regions ['gstart'] + non_unique_regions ['qstart'] #if false
                                            )
    #end coordinates
    non_unique_regions['qend'] = non_unique_regions['gstart'] + non_unique_regions['qend'] 
    non_unique_regions = non_unique_regions.drop(columns=['gstart'])
    non_unique_regions = non_unique_regions.rename(columns={'qseqid': 'Chromosome', 'qstart': 'Start', 'qend': 'End'})
    non_unique_regions['Chromosome'] = 'chr' + non_unique_regions['Chromosome']

    print(f"{' DONE ':=^60}")

    return non_unique_regions

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
    extract_seqs(gsh[0], args.f, 'output/safeharbors.fasta')

    #RM regions
    rm_df = process_rm(args.rm)

    #safe harbors with RM data removed
    gsh_rm = add_regions_to_avoid(gsh[1], [rm_df], args.chro, args.len, 'output/safeharbors_rm.tsv')

    #safe harbors sequences with RM data removed
    extract_seqs(gsh_rm, args.f, 'output/safeharbors_rm.fasta')

    #Blast RM safe harbors
    run_blast('output/safeharbors_rm.fasta', '1e-3', 'output/blast_gsh_rm.txt')

    #Find non-unique regions
    non_unique_df = find_non_unique_seqs('output/blast_gsh_rm.txt')

    #safe harbors with RM data and non-unique regions added to regions to avoid
    gsh_rm_nu = add_regions_to_avoid(gsh[1], [rm_df, non_unique_df], args.chro, args.len, 'output/safeharbors_rm_unique.tsv')

    #get sequences of unique, high information safeharbors
    extract_seqs(gsh_rm_nu, args.f, 'output/safeharbors_rm_unique.fasta')

if __name__ == "__main__":
	main()
