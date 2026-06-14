from Bio import Blast
from Bio import SeqIO
from Bio.Blast import NCBIXML
from io import BytesIO
import time
import sys

Blast.tool = 'biopython'
Blast.email = sys.argv[2]

for record in SeqIO.parse(sys.argv[1], 'fasta'):
    result_stream = Blast.qblast('blastn', 'nt', record.seq, hitlist_size=50, entrez_query='Danio rerio[Organism]')
    
    xml_data = result_stream.read()
    result_stream.close()

    blast_records = NCBIXML.parse(BytesIO(xml_data))
    
    filename = f'{record.id}_blastn.txt'
    
    with open(filename, 'w') as out_stream:
        out_stream.write(f'{"="*60}\n')
        out_stream.write(f'Query: {record.description}\n')
        out_stream.write(f'{"="*60}\n')
    
        for blast_record in blast_records:
            for i, alignment in enumerate(blast_record.alignments, 1):
                hsp = alignment.hsps[0]
                title = alignment.title
                identity = f'{hsp.identities}/{hsp.align_length} ({100*hsp.identities//hsp.align_length}%)'
                out_stream.write(f'\nHit {i}:\n')
                out_stream.write(f'  Title:    {title}\n')
                out_stream.write(f'  E-value:  {hsp.expect:.2e}\n')
                out_stream.write(f'  Score:    {hsp.score}\n')
                out_stream.write(f'  Identity: {identity}\n')
                out_stream.write(f'  Alignment Length: {hsp.align_length}\n')
                start = min(hsp.sbjct_start, hsp.sbjct_end)
                end = max(hsp.sbjct_start, hsp.sbjct_end)
                out_stream.write(f'  Target Start Coordinate: {start}\n')
                out_stream.write(f'  Target End Coordinate: {end}\n')
                out_stream.write(f'  Strand: {hsp.strand}\n')
                out_stream.write(f'\n  -----ALIGNMENT-----\n')
                out_stream.write(f'  Query:  {hsp.query}\n')
                out_stream.write(f'  Match:  {hsp.match}\n')
                out_stream.write(f'  Sbjct:  {hsp.sbjct}\n')

    print(f'Saved: {filename} | Query: {record.description}')
    time.sleep(10)