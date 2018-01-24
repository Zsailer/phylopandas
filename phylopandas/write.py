# Import pandas
import pandas as pd

# Import Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Alphabet


def _write(dataframe, filename=None, schema='fasta', sequence_col='sequence',
           id_col='id', id_only=False, alphabet=None, **kwargs):
    """Write a PhyloPandas DataFrame to a sequence file.
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

    seq_records = []
    # Iterate through dataframe rows
    for i, row in dataframe.iterrows():
        seq = Seq(row[sequence_col], alphabet)
        # Create a SeqRecord and append to list.
        if id_only:
            record = SeqRecord(seq, id=row[id_col], name='', description='')
        else:
            record = SeqRecord(seq, id=row[id_col], name=row['name'],
                               description=row['description'])
        seq_records.append(record)

    # Write to disk or return string
    if filename is not None:
        SeqIO.write(seq_records, filename, format=schema, **kwargs)
    else:
        return "".join([s.format(schema) for s in seq_records])


class BaseWriterAccessor(object):
    """Base Accessor for Pandas Writing attributes.

    To create a new writer.
    """
    # Defaults
    schema = 'fasta'
    id_only = False

    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    def __call__(self, *args, **kwargs):
        """""".format(self.__doc__)
        return _write(self._obj, *args, **kwargs)


@pd.api.extensions.register_dataframe_accessor('to_fasta')
class FastaWriter(BaseWriterAccessor):
    """Write to fasta format.

    Parameters
    ----------
    filename : str
        File to write fasta string to. If no filename is given, a fasta string
        will be returned.
    sequence_col : str (default='sequence')
        Sequence column name in DataFrame.
    id_col : str (default='id')
        ID column name in DataFrame
    id_only : bool (default=False)
        If True, use only the ID column to label sequences in fasta.
    """
    schema = 'fasta'


@pd.api.extensions.register_dataframe_accessor('to_phylip')
class PhylipWriter(BaseWriterAccessor):
    """Phylip Writer"""
    schema = 'phylip'


@pd.api.extensions.register_dataframe_accessor('to_clustal')
class ClustalWriter(BaseWriterAccessor):
    """Clustal Writer"""
    schema = 'clustal'


@pd.api.extensions.register_dataframe_accessor('to_embl')
class EmblWriter(BaseWriterAccessor):
    """Embl Writer"""
    schema = 'embl'


@pd.api.extensions.register_dataframe_accessor('to_nexus')
class NexusWriter(BaseWriterAccessor):
    """Nexus Writer"""
    schema = 'nexus'


@pd.api.extensions.register_dataframe_accessor('to_swiss')
class SwissWriter(BaseWriterAccessor):
    """Swiss Writer"""
    schema = 'swiss'
