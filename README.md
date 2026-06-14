# Genomic-Safeharbors-Danio-Rerio-Zebrafish

### Overview
The method for finding genomic safeharbors computationally for *Danio rerio* (Zebrafish) follows the workflow outline in [Discovery and validation of human genomic safe harbor sites for gene and cell therapies](https://pmc.ncbi.nlm.nih.gov/articles/PMC9017210/).

### Dataset 

| File Name | Datatype | Assembly | Download Location | Coordinate Type | Notes |
| --------- | -------- | -------- |------------------ | --------------- | ----- |
| `danRer11_chromL.txt` | Lengths of the chromsomes | v11 | [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/) | NA |
| `danRer11_gene.gtf`   | Gene coordinates | v11 | [Ensembl Database](http://ensembl.org/) | GTF 
| `danRer11_seq.fa`     | Primary assembly | v11 | [Ensembl Database](http://ensembl.org/) |  | 
| `danRer11_enh.bed`    | Enhancer coordinates | v10 → v11 | [EnhancerAtlas2.0](https://ngdc.cncb.ac.cn/databasecommons/database/id/7011) | BED | Convert to v11 using [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
| `danRer11_lnc.bed`    | lncRNA coordinates | v7 → v11  | [ZFLS](https://old.biochen.org/zflnc/) | BED | Convert to BED format then use [Liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert to v11 assembly | 
| `danRer11_t.gtf`      | tRNA coordinates | v11 | [GtRNADB](https://gtrnadb.org/) | GTF | Raw download may have malformed columns |
| `danRer11_mi.bed`     | miRNA coordiantes | v11 | [MirGeneDB 3.0](https://mirgenedb.org/) | BED | 
| `danRer11_onco.txt`   | Gene names of oncogenes | v11 | [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v104/cancergenecensus) & [Human Zebrafish Gene Orthologs - ZFIN](https://zfin.org/downloads/human_orthos.txt) |  | Found ortholog genes between human and zebrafish then extracted those genes names from the COSMIC Cancer Gene list |
| `danRer11_gap.txt`    | Gap coordinates | v11 |  [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables) | BED | 
| `danRer11_cent.gtf`   | Centromere coordinates | v11 |  | GTF | Generated using [Quartet CentroMiner](https://github.com/aaranyue/quarTeT) and top scoring centromeres extracted manually |
| `danRer11_rm.bed` | RepeatMasker Coordinates | v11 | [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?db=danRer11&hgta_group=varRep&hgta_track=rmsk&hgta_table=rmsk) |  BED | 

### Requirements
`pandas
pybedtools
pyranges==0.0.129
pyfaidx
numpy`

### Usage
`python gsh_python.py \
  -chro data/danRer11_chromL.txt \
  -genes data/danRer11_gene.gtf \
  -onco data/danRer11_onco.txt \
  -enh data/danRer11_enh.bed \
  -cent data/danRer11_cent.gtf \
  -gap data/danRer11_gap.txt \
  -lnc data/danRer11_lnc.bed \
  -mi data/danRer11_mi.bed \
  -t data/danRer11_t.gtf \
  -rm data/danRer11_rm.bed \
  -f data/danRer11_seq.fa`
  
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
| safeharbors.tsv | Genomic coordinates of safeharbors with columns `Chromosome`, `Start` and `End` |
| safeharbors_seqs.fasta | Sequences of the safeharbors |
| safeharbors_rm.tsv | Genomic coordinates of safeharbors with RepeatMasker coordinates subtracted |
| safeharbors_rm_seqs.fasta | Sequences of the safeharbors with RepeatMasker coordinates subtracted |
