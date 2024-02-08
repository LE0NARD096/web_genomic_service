import os
from django.core.management.base import BaseCommand, CommandError
from website.models import GeneProtein, Genome, AnnotationGenome, AnnotationProtein
from Bio import SeqIO
import re
from django.utils import timezone
import warnings
from django.db import transaction

class Command(BaseCommand):
    help = 'Command to save files in the database'

    def add_arguments(self, parser):
        # Define command-line arguments
        parser.add_argument('folder_path', type=str, help='Folder path containing files to be saved in the database')
        parser.add_argument('--annotated', action='store_true' , help='The sequences are all annotated or not')
    

    def handle(self, *args, **options):
        folder_path = options['folder_path']
        annotated = options['annotated']
    
        files_fasta=[]
        for file_fasta in os.listdir(folder_path):
                if file_fasta.lower().endswith('.fa'):
                    files_fasta.append(file_fasta)

        
        pattern = re.compile(r'([a-zA-Z0-9_]+)_(cds|pep|[a-z0-9]*)\.fa$')

        for species_and_type in files_fasta:
            match = pattern.match(species_and_type)
            print(species_and_type)

            if match.group(1) and (match.group(2) == 'cds' or match.group(2) == 'pep'):
                sequence_proteins = []
                annotation_proteins = []
                counter_for_bulk = 0
                base = 900

                with transaction.atomic():
                    c = list(SeqIO.parse(os.path.join(folder_path, species_and_type), "fasta"))
                    for record in SeqIO.parse(os.path.join(folder_path, species_and_type), 'fasta'):
                        
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
                                                            defaults={
                                                                'upload_time':timezone.now()
                                                                    })

                            if created and genome.sequence is None:
                                print('error 1')
                                genome.save()
                            
                            if annotated:
                                validated_status = True
                        
                            geneprotein_instance = GeneProtein(
                                                                                accession_number=a_n,
                                                                                sequence=sequence,
                                                                                type=type,
                                                                                start=start,
                                                                                end=end,
                                                                                genome=genome,
                                                                                is_validated = validated_status,
                                                                                upload_time=timezone.now()
                                                                                )
                            
                            sequence_proteins.append(geneprotein_instance)
                            
                            if annotated:
                                gene=description[9]
                                transcript=a_n
                                gene_biotype=description[12]
                                gene_symbol=description[15]
                                transcript_biotype=description[14]

                                if type == "cds":
                                    description_protein=' '.join(description[18:21])
                                else:
                                    description_protein=' '.join(description[19:21])

                                annotation_protein = AnnotationProtein(
                                                                        gene=gene,
                                                                        transcript=transcript,
                                                                        gene_biotype=gene_biotype,
                                                                        transcript_biotype=transcript_biotype,
                                                                        gene_symbol=gene_symbol,
                                                                        geneprotein=geneprotein_instance,
                                                                        description=description_protein,
                                                                        is_annotated=True,
                                                                        annotation_time=timezone.now()
                                                                    )

                                annotation_proteins.append(annotation_protein)

                            counter_for_bulk += 1

                            if counter_for_bulk == base or counter_for_bulk == len(c):
                                    
                                base += 900

                                if annotated:
                                    GeneProtein.objects.bulk_create(sequence_proteins,batch_size=900)
                                    AnnotationProtein.objects.bulk_create(annotation_proteins,batch_size=900)
                                                
                                else:
                                    GeneProtein.objects.bulk_create(sequence_proteins, batch_size=900)
                                
                                sequence_proteins = []
                                annotation_proteins = []
                                    

            elif match.group(1):
                try:
                    genome = SeqIO.read(os.path.join(folder_path, species_and_type), 'fasta')
                except:
                    raise CommandError(f'Genome fasta file {species_and_type} expects exactly one fasta sequence')
  

                sequence = genome.seq
                description = genome.description.split(':')
                chromosome = description[2]

                if not sequence:
                    raise CommandError(f'The genome fasta file {species_and_type} does not contain a sequence')

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
                    genome.save()
                
                else:
                    warnings.warn(f'The genome in fasta file {species_and_type} is already in the database', stacklevel=2)
                    continue
                               
                if annotated:
                    species = species = re.split(r'[_.]', species_and_type)[:-1]
                    species = " ".join(species)
                    annotations = AnnotationGenome.objects.create(species=species,
                                                                    genome=genome,
                                                                    is_annotated=True,
                                                                    annotation_time = timezone.now())
                    annotations.save()
            
            else:
                raise CommandError(f"The title of the fasta file {species_and_type} doesn't respect your_bacteria_species_cds or pep.fa")

            
        
