import os
from django.core.management.base import BaseCommand, CommandError
from website.models import GeneProtein, Genome, AnnotationGenome, AnnotationProtein
from Bio import SeqIO
import re
from django.utils import timezone
import warnings
from django.db import transaction
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


class Command(BaseCommand):
    help = 'Command to save files in the database'

    def add_arguments(self, parser):

        parser.add_argument('folder_path', type=str, help='Folder path containing files to be saved in the database')
        parser.add_argument('query_per_sequence', type=int, help='Number of sequences saved in the database per query', choices=range(1, 999))
        parser.add_argument('--annotated', action='store_true' , help='The sequences are all annotated or not')
    
    def handle(self, *args, **options):
        folder_path = options['folder_path']
        annotated = options['annotated']
        query_per_sequence = options['query_per_sequence']
    
        files_fasta=[]
        for file_fasta in os.listdir(folder_path):
                if file_fasta.lower().endswith('.fa'):
                    files_fasta.append(file_fasta)
        
    
        
        pattern = re.compile(r'([a-zA-Z0-9_]+)_(cds|pep|[a-z0-9]*)\.fa$')

        for file in files_fasta:
            match = pattern.match(file)
            print(file)

            if match.group(1) and (match.group(2) == 'cds' or match.group(2) == 'pep'):
                sequence_proteins = []
                annotation_proteins = []
                counter_for_bulk = 0
                base = query_per_sequence

                sum_sequences_file = list(SeqIO.parse(os.path.join(folder_path, file), "fasta"))
                record_iter = SeqIO.parse(open(os.path.join(folder_path, file)), "fasta")
                
                with transaction.atomic():
                    for i, batch in enumerate(batch_iterator(record_iter, query_per_sequence)):
                       for record in batch:   
                            a_n = record.id
                            sequence = record.seq

                            if not sequence:
                                warnings.warn(f"The protein {a_n} wasn't inserted into the database. Reason: No sequence", stacklevel=2)
                                continue
                

                            type = record.description.split()[1].lower()
                            description = re.split(r'[: ]', record.description)
                            start = description[5]
                            end = description[6]
                            chromosome = description[3]
                            validated_status = False


                            genome, created = Genome.objects.get_or_create(chromosome=chromosome,
                                                            defaults={'upload_time':timezone.now()})

                            if created and genome.sequence is None:
                                genome.save()
                            
                            if annotated:
                                validated_status = True
                        
                            geneprotein_instance = GeneProtein(accession_number=a_n,
                                                               sequence=sequence,
                                                               type=type,
                                                               start=start,
                                                               end=end,
                                                               genome=genome,
                                                               is_validated = validated_status)
                                                                                
                            
                            sequence_proteins.append(geneprotein_instance)
                            
                            if annotated:
                                gene=description[9]
                                gene_symbol_exist = re.search(r'gene_symbol:(\S+)',record.description)

                                if gene_symbol_exist:
                                    if type == "cds":
                                        transcript = "NA"
                                        gene_biotype = description[11]
                                        transcript_biotype = description[13]
                                        gene_symbol=description[15]
                                        description_protein=' '.join(description[17:])
                                    else:
                                        transcript = a_n
                                        gene_biotype = description[13]
                                        transcript_biotype = description[15]
                                        gene_symbol=description[17]
                                        description_protein=' '.join(description[19:])
                            
                                else:
                                    gene_symbol="unknown"

                                    if type == "cds":
                                        transcript = "NA"
                                        gene_biotype = description[11]
                                        transcript_biotype = description[13]
                                        description_protein=' '.join(description[15:])
                                    else:
                                        transcript = a_n
                                        gene_biotype = description[13]
                                        transcript_biotype = description[15]
                                        description_protein=' '.join(description[17:])


                                annotation_protein = AnnotationProtein(
                                                                        gene=gene,
                                                                        transcript=transcript,
                                                                        gene_biotype=gene_biotype,
                                                                        transcript_biotype=transcript_biotype,
                                                                        gene_symbol=gene_symbol,
                                                                        geneprotein=geneprotein_instance,
                                                                        description=description_protein,
                                                                        is_annotated=True
                                                                    )

                                annotation_proteins.append(annotation_protein)

                            counter_for_bulk += 1

                            if counter_for_bulk == base or counter_for_bulk == len(sum_sequences_file):
                                    
                                base += query_per_sequence

                                
                                if annotated:

                                    for i in range(len(sequence_proteins)):
                                        sequence_proteins[i].upload_time = timezone.now()
                                        annotation_proteins[i].annotation_time = timezone.now()

                                    GeneProtein.objects.bulk_create(sequence_proteins)
                                    AnnotationProtein.objects.bulk_create(annotation_proteins)
                                                
                                else:
                                    for i in range(len(sequence_proteins)):
                                        sequence_proteins[i].upload_time = timezone.now()

                                    GeneProtein.objects.bulk_create(sequence_proteins)
                                
                                sequence_proteins.clear()
                                annotation_proteins.clear()
                                    

            elif match.group(1):
                try:
                    genome = SeqIO.read(os.path.join(folder_path, file), 'fasta')
                except:
                    raise CommandError(f'Genome fasta file {file} expects exactly one fasta sequence')
  

                sequence = genome.seq
                description = genome.description.split(':')
                chromosome = description[2]

                if not sequence:
                    raise CommandError(f'The genome fasta file {file} does not contain a sequence')

                start = description[4]
                end = description[5]
                    
                genome, created = Genome.objects.get_or_create(chromosome=chromosome,
                                                    defaults={
                                                        'sequence':sequence,
                                                        'start':start,
                                                        'end':end,
                                                        'upload_time':timezone.now()
                                                            })
                                                
                # We populate the "artificial" genome created by the proteins
                if not created and genome.sequence is None:
                    genome.sequence=sequence
                    genome.start=start
                    genome.end=end
                    genome.is_validated=True
                    genome.save()
                
                elif created:
                    genome.is_validated=True
                    genome.save()
                
                else:
                    warnings.warn(f'The genome in fasta file {file} is already in the database', stacklevel=2)
                    continue
                               
                if annotated:
                    species = species = re.split(r'[_.]', file)[:-1]
                    species = " ".join(species)
                    annotations = AnnotationGenome.objects.create(species=species,
                                                                    genome=genome,
                                                                    is_annotated=True,
                                                                    annotation_time = timezone.now())
                    annotations.save()
            
            else:
                raise CommandError(f"The title of the fasta file {file} doesn't respect your_bacteria_species_cds or pep.fa")

            
        
