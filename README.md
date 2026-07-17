# Genomic-Safeharbors-Danio-Rerio-Zebrafish

### Overview
The method for finding genomic safeharbors (GSH) computationally for *Danio rerio* (Zebrafish) follows the workflow outline in [Discovery and validation of human genomic safe harbor sites for gene and cell therapies](https://pmc.ncbi.nlm.nih.gov/articles/PMC9017210/). Two additional steps have been added to this pipeline. Once the GSH canadiates have been computationally generated, low complexity regions and non-unique sequences are excluded.
* **Low complexity regions** are excluded by adding [RepeatMasker](https://www.repeatmasker.org/) regions to the regions to avoid. RepeatMasker data for *Danio rerio* was downloaded from [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?db=danRer11&hgta_group=varRep&hgta_track=rmsk&hgta_table=rmsk) (see Dataset table below).
* **Non-unique sequences** were excluded by running a local BLAST on all GSH canadiate sequences and adding the regions of said canadiates that returned a BLAST hit (not including self hits) to the regions to avoid. This step is meant to narrow down the GSH canadiates to only regions with unique sequences, which allows for precise targeting. 
### Dataset 
All cleaned and processed data is avaliable to download from here: [danioRer_data](https://iastate.box.com/s/njar9ckjgxg75in7k08fyxxaep3geo2q) 

| File Name | Datatype | Assembly | Download Location | Coordinate Type | Notes |
| --------- | -------- | -------- |------------------ | --------------- | ----- |
| `danRer11_chromL.txt` | Lengths of the chromsomes | v11 | [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/) | NA |
| `danRer11_gene.gtf`   | Gene coordinates | v11 | [Ensembl Database](https://ftp.ensembl.org/pub/release-116/gtf/danio_rerio/) | GTF 
| `danRer11_seq.fa`     | Whole genome sequence | v11 | [Ensembl Database](https://ftp.ensembl.org/pub/release-116/fasta/danio_rerio/dna/) |  | 
| `danRer11_enh.bed`    | Enhancer coordinates | v10 → v11 | [EnhancerAtlas2.0](https://ngdc.cncb.ac.cn/databasecommons/database/id/7011) | BED | Converted to v11 using [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
| `danRer11_lnc.bed`    | lncRNA coordinates | v7 → v11  | [ZFLS](https://old.biochen.org/zflnc/) | BED | Converted to BED format then used [Liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert to v11 assembly | 
| `danRer11_t.gtf`      | tRNA coordinates | v11 | [GtRNADB](https://gtrnadb.org/) | GTF | Raw download may have malformed columns |
| `danRer11_mi.bed`     | miRNA coordiantes | v11 | [MirGeneDB 3.0](https://mirgenedb.org/) | BED | 
| `danRer11_onco.txt`   | Gene names of oncogenes | v11 | [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v104/cancergenecensus) & [Human Zebrafish Gene Orthologs - ZFIN](https://zfin.org/downloads/human_orthos.txt) |  | Found ortholog genes between human and zebrafish then extracted those genes names from the COSMIC Cancer Gene list |
| `danRer11_gap.txt`    | Gap coordinates | v11 |  [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=4092619007_GihEf7a5LhdatZ3HK8bmgfxMGSkR&db=danRer11&hgta_group=map&hgta_track=gap&hgta_table=0&hgta_regionType=range&position=chr5%3A44%2C806%2C795-48%2C499%2C271&hgta_outputType=primaryTable&hgta_outFileName=) | BED | 
| `danRer11_cent.gtf`   | Centromere coordinates | v11 |  | GTF | Generated using [Quartet CentroMiner](https://github.com/aaranyue/quarTeT) and top scoring centromeres extracted manually |
| `danRer11_rm.bed` | RepeatMasker Coordinates | v11 | [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?db=danRer11&hgta_group=varRep&hgta_track=rmsk&hgta_table=rmsk) |  BED | 

### Requirements
`pandas`
`pybedtools`
`pyranges==0.0.129`
`pyfaidx`
`numpy`

### Local BLAST
A local *Danio rerio* BLAST database is required for this code to run. 
1) Follow instructions to install command line BLAST from the [NIH](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
2) Generate database from Ensembl primary assembly file (`danRer11_seq.fa` see Dataset table above):
   `makeblastdb -in danRer_data/danRer11_seq.fa -input_type fasta -dbtype nucl -parse_seqids -out danRer11_db/danRer11`
   
### Usage
* The folder `danRer_data` (downloaded from the box [danioRer_data](https://iastate.box.com/s/njar9ckjgxg75in7k08fyxxaep3geo2q)) needs to be in the same folder as gsh_python.py
> `python gsh_python.py -chro danRer_data/danRer11_chromL.txt -genes danRer_data/danRer11_gene.gtf -onco danRer_data/danRer11_onco.txt -enh danRer_data/danRer11_enh.bed -cent danRer_data/danRer11_cent.gtf -gap danRer_data/danRer11_gap.txt -lnc danRer_data/danRer11_lnc.bed -mi danRer_data/danRer11_mi.bed -t danRer_data/danRer11_t.gtf -rm danRer_data/danRer11_rm.bed -f danRer_data/danRer11_seq.fa -len 500`
  
### Flank Distances
| Feature | Distance (bp) |
| ------- | ------------- |
| Gene | 5,000 |
| Oncogene | 300,000 |
| miRNA | 300,000 |
| tRNA | 150,000 |
| lncRNA | 150,000 |
| Enhancer | 20,000 |
| Centromere | 300,000 |
| Gap | 300,000 |

### Outputs
| Filename | Description |
| -------- | ----------- |
| safeharbors.tsv | Genomic coordinates of safeharbors with columns `Chromosome`, `Start`, `End`, `Length` and `LP` |
| safeharbors_seqs.fasta | Sequences of the safeharbors |
| safeharbors_rm.tsv | Genomic coordinates of safeharbors with RepeatMasker coordinates added to regions to avoid (columns `Chromosome`, `Start`, `End`, `Length` and `LP`) |
| safeharbors_rm_seqs.fasta | Sequences of the safeharbors with RepeatMasker coordinates added to regions to avoid |
| safeharbors_rm_unique.tsv | Genomic coordinates of safeharbors with RepeatMasker and non-unique region coordinates added to regions to avoid (columns `Chromosome`, `Start`, `End`, `Length` and `LP`) |
| safeharbors_rm_unique.fasta | Sequences of the safeharbors with RepeatMasker and non-unique region coordinates added to regions to avoid |

### Citations
[Aznauryan et al. (2022). "Discovery and validation of human genomic safe harbor sites for gene and cell therapies." *Cell Reports Methods*.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9017210/)
