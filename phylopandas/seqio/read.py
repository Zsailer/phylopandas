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
import pandas as pd
from ..utils import get_random_id


def _read_doc_template(schema):
    s = """Read a {} file.

    Construct a PhyloPandas DataFrame with columns:
        - name
        - id
        - description
        - sequence

    Parameters
    ----------
    filename : str
        File name of {} file.

    seq_label : str (default='sequence')
        Sequence column name in DataFrame.
    """.format(schema, schema, schema)
    return s


def _read(
    filename,
    schema,
    seq_label='sequence',
    alphabet=None,
    use_uids=True,
    **kwargs):
    """Use BioPython's sequence parsing module to convert any file format to
    a Pandas DataFrame.

    The resulting DataFrame has the following columns:
        - name
        - id
        - description
        - sequence
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
    data = {
        'id': [],
        seq_label: [],
        'description': [],
        'label': []
    }
    if use_uids:
        data['uid'] = []

    # Parse Fasta file.
    for i, s in enumerate(SeqIO.parse(filename, format=schema, **kwargs)):
        data['id'].append(s.id)
        data[seq_label].append(str(s.seq))
        data['description'].append(s.description)
        data['label'].append(s.name)

        if use_uids:
            data['uid'].append(get_random_id(10))

    # Port to DataFrame.
    return pd.DataFrame(data)


def _read_method(schema):
    """Add a write method for named schema to a class.
    """
    def func(
        self,
        filename,
        seq_label='sequence',
        alphabet=None,
        combine_on='uid',
        use_uids=True,
        **kwargs):
        # Use generic write class to write data.
        df0 = self._data
        df1 = _read(
            filename=filename,
            schema=schema,
            seq_label=seq_label,
            alphabet=alphabet,
            use_uids=use_uids,
            **kwargs
        )
        return df0.phylo.combine(df1, on=combine_on)

    # Update docs
    func.__doc__ = _read_doc_template(schema)
    return func


def _read_function(schema):
    """Add a write method for named schema to a class.
    """
    def func(
        filename,
        seq_label='sequence',
        alphabet=None,
        use_uids=True,
        **kwargs):
        # Use generic write class to write data.
        return _read(
            filename=filename,
            schema=schema,
            seq_label=seq_label,
            alphabet=alphabet,
            use_uids=use_uids,
            **kwargs
        )
    # Update docs
    func.__doc__ = _read_doc_template(schema)
    return func


# Various read functions to various formats.
read_fasta = _read_function('fasta')
read_phylip = _read_function('phylip')
read_clustal = _read_function('clustal')
read_embl = _read_function('embl')
read_nexus = _read_function('nexus')
read_swiss = _read_function('swiss')
read_fastq = _read_function('fastq')


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
    return pd.DataFrame(data)
