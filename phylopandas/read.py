__doc__ = """
Functions for reading sequence files into pandas DataFrame.
"""

# Imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import Bio.Alphabet

# Import Phylopandas DataFrame
from pandas import DataFrame


def _read(filename, schema, seq_label='sequence', alphabet=None, **kwargs):
    """Use BioPython's sequence parsing module to convert any file format to
    a Pandas DataFrame.
    """
    # Check Alphabet if given
    if alphabet is None:
        alphabet = Bio.Alphabet.Alphabet()
    elif alphabet in ['dna', 'rna', 'protein', 'nucleotide']:
        alphabet = getattr(Bio.Alphabet, 'generic_{}'.format(alphabet))
    else:
        raise Exception(
            "The alphabet is not recognized. Must be 'dna', 'rna', "
            "'nucleotide', or 'protein'.")

    kwargs.update(alphabet=alphabet)

    # Prepare DataFrame fields.
    data = {'id': [], seq_label: [], 'description': [], 'name': []}

    # Parse Fasta file.
    for i, s in enumerate(SeqIO.parse(filename, format=schema, **kwargs)):
        data['id'].append(s.id)
        data[seq_label].append(str(s.seq))
        data['description'].append(s.description)
        data['name'].append(s.name)

    # Port to DataFrame.
    return DataFrame(data)


def read_fasta(filename, **kwargs):
    """Read fasta format."""
    return _read(filename, schema='fasta', **kwargs)


def read_phylip(filename, **kwargs):
    """Read phylip format."""
    return _read(filename, schema='phylip', **kwargs)


def read_clustal(filename, **kwargs):
    """Read clustal format."""
    return _read(filename, schema='clustal', **kwargs)


def read_embl(filename, **kwargs):
    """Read the EMBL flat file format."""
    return _read(filename, schema='embl', **kwargs)


def read_nexus(filename, **kwargs):
    """Read the EMBL flat file format."""
    return _read(filename, schema='nexus', **kwargs)


def read_swiss(filename, **kwargs):
    """Read Swiss-Prot aka UniProt format."""
    return _read(filename, schema='nexus', **kwargs)


def read_fastq(filename, **kwargs):
    """Read FASTQ format."""
    return _read(filename, schema='fastq', **kwargs)


def read_blast_xml(filename, **kwargs):
    """Read BLAST XML format."""
    # Read file.
    with open(filename, 'r') as f:
        blast_record = NCBIXML.read(f)

    # Prepare DataFrame fields.
    data = {'accession': [],
            'hit_def': [],
            'hit_id': [],
            'title': [],
            'length': [],
            'e_value': [],
            'sequence': []}

    # Get alignments from blast result.
    for i, s in enumerate(blast_record.alignments):
        data['accession'] = s.accession
        data['hit_def'] = s.hit_def
        data['hit_id'] = s.hit_id
        data['title'] = s.title
        data['length'] = s.length
        data['e_value'] = s.hsps[0].expect
        data['sequence'] = s.hsps[0].sbjct

    # Port to DataFrame.
    return DataFrame(data)
