
from django.core.management.base import BaseCommand, CommandError
from website.models import GeneProtein, Genome, AnnotationGenome, AnnotationProtein, Profile
from django.utils import timezone
from django.db import transaction
from django.contrib.auth.hashers import make_password

from Bio import SeqIO
import re
import warnings
import os


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
        parser.add_argument('--annotated', action='store_true', help='The sequences are all annotated or not')
        parser.add_argument('--populatewithusers', action='store_true', help='Populate the site with test users')

    def handle(self, *args, **options):
        folder_path = options['folder_path']
        annotated = options['annotated']
        populate_with_users = options['populatewithusers']

        if populate_with_users:
            created_user = Profile.objects.get_or_create(email='user@gmail.com',
                                                         defaults={'username': 'user',
                                                                   'role': 'user',
                                                                   'password': make_password('test'),
                                                                   'is_approved': True})
            if created_user[1]:
                print('User created')

            created_annotator = Profile.objects.get_or_create(email='annotator@gmail.com',
                                                              defaults={'username': 'annotator',
                                                                        'role': 'annotator',
                                                                        'password': make_password('test'),
                                                                        'is_approved': True})
            if created_annotator[1]:
                print('Annotator created')

            created_validator = Profile.objects.get_or_create(email='validator@gmail.com',
                                                              defaults={'username': 'validator',
                                                                        'role': 'validator',
                                                                        'password': make_password('test'),
                                                                        'is_approved': True})

            if created_validator[1]:
                print('Validator created')

            if not Profile.objects.filter(email='admin@gmail.com').exists():
                admin = Profile.objects.create_superuser('admin', 'admin@gmail.com', 'test')
                admin.role = 'admin'
                admin.is_approved = True
                admin.save()
                print('Admin created')

        files_fasta = []
        for file_fasta in os.listdir(folder_path):
            if file_fasta.lower().endswith('.fa'):
                files_fasta.append(file_fasta)

        pattern = re.compile(r'([a-zA-Z0-9_]+)_(cds|pep|[a-z0-9]*)\.fa$')

        for file in files_fasta:
            match = pattern.match(file)
            print(file)

            if match.group(1) and (match.group(2) == 'cds' or match.group(2) == 'pep'):
                query_per_sequence = 800
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
                                                                           defaults={'upload_time': timezone.now()})

                            pattern2 = re.compile(r'description:(.*?)\b')
                            match = pattern2.search(record.description)

                            if annotated and match:
                                validated_status = True

                            geneprotein_instance, created = GeneProtein.objects.get_or_create(accession_number=a_n,
                                                                                              type=type,
                                                                                              defaults={'start': start,
                                                                                                        'sequence': sequence,
                                                                                                        'end': end,
                                                                                                        'genome': genome,
                                                                                                        'is_validated': validated_status,
                                                                                                        'upload_time': timezone.now()})

                            if not created:
                                continue

                            if annotated and match:
                                gene = description[9]
                                gene_symbol_exist = re.search(r'gene_symbol:(\S+)', record.description)

                                if gene_symbol_exist:
                                    if type == "cds":
                                        transcript = "NA"
                                        gene_biotype = description[11]
                                        transcript_biotype = description[13]
                                        gene_symbol = description[15]
                                        description_protein = ' '.join(description[17:])
                                    else:
                                        transcript = a_n
                                        gene_biotype = description[13]
                                        transcript_biotype = description[15]
                                        gene_symbol = description[17]
                                        description_protein = ' '.join(description[19:])

                                else:
                                    gene_symbol = "unknown"

                                    if type == "cds":
                                        transcript = "NA"
                                        gene_biotype = description[11]
                                        transcript_biotype = description[13]
                                        description_protein = ' '.join(description[15:])
                                    else:
                                        transcript = a_n
                                        gene_biotype = description[13]
                                        transcript_biotype = description[15]
                                        description_protein = ' '.join(description[17:])

                                    AnnotationProtein.objects.create(gene=gene,
                                                                     transcript=transcript,
                                                                     gene_biotype=gene_biotype,
                                                                     transcript_biotype=transcript_biotype,
                                                                     gene_symbol=gene_symbol,
                                                                     geneprotein=geneprotein_instance,
                                                                     description=description_protein,
                                                                     annotation_time=timezone.now(),
                                                                     is_annotated=True)

            elif match.group(1):
                try:
                    genome = SeqIO.read(os.path.join(folder_path, file), 'fasta')
                except Exception as e:
                    raise CommandError(f'Error in {file}: {e}')

                sequence = genome.seq
                description = genome.description.split(':')
                chromosome = description[2]

                if not sequence:
                    raise CommandError(f'The genome fasta file {file} does not contain a sequence')

                start = description[4]
                end = description[5]

                genome, created = Genome.objects.get_or_create(chromosome=chromosome,
                                                               defaults={'sequence': sequence,
                                                                         'start': start,
                                                                         'end': end,
                                                                         'upload_time': timezone.now()})

                pattern2 = re.compile(r'new_(.*?)\b')
                match = pattern2.search(file)

                if not created and genome.sequence is None and match is None:
                    genome.sequence = sequence
                    genome.start = start
                    genome.end = end
                    genome.is_validated = True
                    genome.save()
                elif created and match is None:
                    genome.is_validated = True
                    genome.save()
                elif match:
                    genome.is_validated = False
                    genome.save()

                else:
                    warnings.warn(f'The genome in fasta file {file} is already in the database', stacklevel=2)
                    continue
                if annotated and match is None:
                    species = re.split(r'[_.]', file)[:-1]
                    species = " ".join(species)
                    AnnotationGenome.objects.create(species=species,
                                                    genome=genome,
                                                    is_annotated=True,
                                                    annotation_time=timezone.now())
            else:
                raise CommandError(f"The title of the fasta file {file} doesn't respect your_bacteria_species_cds or pep.fa")
