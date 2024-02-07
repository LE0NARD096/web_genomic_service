import os
from django.core.management.base import BaseCommand, CommandError
from website.models import GeneProtein, Genome, AnnotationGenome, AnnotationProtein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO
import re
from django.utils import timezone
import warnings

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

        
        pattern = re.compile(r'^([a-zA-Z0-9_]+)_(cds|pep)?')

        for species_and_type in files_fasta:
            match = pattern.match(species_and_type)
            print(match.group(1))
            
            if match.group(2) == 'cds' or match.group(2) == 'pep':
                 print('protein',species_and_type)
            
            elif match.group(1) and match.group(2) == None:
                print('genome',match.group(1))
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
                    
                genome = Genome.objects.get_or_create(chromosome=chromosome,
                                                    defaults={
                                                        'sequence':sequence,
                                                        'start':start,
                                                        'end':end,
                                                        'upload_time':timezone.now()
                                                            })
                                                
                # We populate the "artificial" genome created by the proteins
                if not genome[1] and genome[0].sequence is None:
                    genome[0].sequence=sequence
                    genome[0].start=start
                    genome[0].end=end
                    genome[0].save()
                
                elif genome[1]:
                    genome[0].save()
                
                else:
                    warnings.warn(f'The genome in fasta file {species_and_type} is already in the database', stacklevel=2)
                    continue
                               
                if annotated:
                    species = species = re.split(r'[_.]', species_and_type)[:-1]
                    species = " ".join(species)
                    annotations = AnnotationGenome.objects.create(species=species,
                                                                    genome=genome[0],
                                                                    is_annotated=True,
                                                                    annotation_time = timezone.now())
                    annotations.save()

            
        
