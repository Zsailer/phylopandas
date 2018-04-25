__doc__ = """
Functions for write sequence data to sequence files.
"""
import pandas as pd

# Import Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Alphabet


def _seqio_doc_template(schema):
    s = """Write to {} format.

    Parameters
    ----------
    filename : str
        File to write {} string to. If no filename is given, a fasta string
        will be returned.
    sequence_col : str (default='sequence')
        Sequence column name in DataFrame.
    """.format(schema, schema)
    return s


def _pandas_series_to_biopython_record(
    series,
    id_col='id',
    id_only=False,
    sequence_col='sequence',
    alphabet=None):
    """
    """
    seq = Seq(series[sequence_col], alphabet)
    # Create a SeqRecord and append to list.
    if id_only:
        record = SeqRecord(seq, id=series[id_col], name='',
                           description='')

    else:
        record = SeqRecord(seq, id=series[id_col], name=series['name'],
                           description=series['description'])
    return record


def _write(
    data,
    filename=None,
    schema='fasta',
    sequence_col='sequence',
    id_col='id',
    id_only=False,
    alphabet=None,
    **kwargs):
    """General write function. Write phylopanda data to biopython format.

    Parameters
    ----------
    filename : str
        File to write string to. If no filename is given, a string
        will be returned.

    sequence_col : str (default='sequence')
        Sequence column name in DataFrame.

    id_col : str (default='id')
        ID column name in DataFrame

    id_only : bool (default=False)
        If True, use only the ID column to label sequences in fasta.
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

    # Build a list of records from a pandas DataFrame
    if type(data) is pd.DataFrame:
        for i, row in data.iterrows():
            record = _pandas_series_to_biopython_record(
                row,
                id_col=id_col,
                id_only=id_only,
                sequence_col=sequence_col,
                alphabet=alphabet
            )
            seq_records.append(record)

    # Build a record from a pandas Series
    elif type(data) is pd.Series:
        record = _pandas_series_to_biopython_record(
            row,
            id_col=id_col,
            id_only=id_only,
            sequence_col=sequence_col,
            alphabet=alphabet
        )
        seq_records.append(record)

    # Write to disk or return string
    if filename is not None:
        SeqIO.write(seq_records, filename, format=schema, **kwargs)

    else:
        return "".join([s.format(schema) for s in seq_records])


def to_fasta(df, filename=None, sequence_col='sequence',
             id_col='id', id_only=False, alphabet=None, **kwargs):
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
    return _write(df, filename=filename, schema='fasta',
                  sequence_col=sequence_col, id_col=id_col, id_only=id_only,
                  alphabet=None, **kwargs)


def to_phylip(df, filename=None, sequence_col='sequence',
             id_col='id', alphabet=None, **kwargs):
    __doc__ = _seqio_doc_template('phylip')
    return _write(df, filename=filename, schema='phylip',
                  sequence_col=sequence_col, id_col=id_col, id_only=True,
                  alphabet=None, **kwargs)


def to_clustal(df, filename=None, sequence_col='sequence',
             id_col='id', alphabet=None, **kwargs):
    __doc__ = _seqio_doc_template('clustal')
    return _write(df, filename=filename, schema='clustal',
                  sequence_col=sequence_col, id_col=id_col, id_only=True,
                  alphabet=None, **kwargs)

def to_embl(df, alphabet, filename=None, sequence_col='sequence',
             id_col='id', **kwargs):
    __doc__ = _seqio_doc_template('embl')
    return _write(df, filename=filename, schema='embl', sequence_col=sequence_col,
                  id_col=id_col, id_only=True, alphabet=alphabet, **kwargs)


def to_nexus(df, alphabet, filename=None, sequence_col='sequence',
             id_col='id', id_only=False, **kwargs):
    __doc__ = _seqio_doc_template('nexus')
    return _write(df, alphabet=alphabet, filename=filename, schema='nexus',
                  sequence_col=sequence_col, id_col='id',
                  id_only=True, **kwargs)


def to_swiss(df, filename=None, sequence_col='sequence',
             id_col='id', id_only=False, alphabet=None, **kwargs):
    __doc__ = _seqio_doc_template('swiss')
    return _write(df, alphabet=alphabet, filename=filename, schema='swiss',
                  sequence_col=sequence_col, id_col='id', id_only=True,
                  **kwargs)


def to_fastq(df, filename=None, sequence_col='sequence',
             id_col='id', id_only=False, alphabet=None, **kwargs):
    __doc__ = _seqio_doc_template('fastq')
    return _write(df, filename=filename, schema='fastq',
                  sequence_col=sequence_col, id_col='id', id_only=True,
                  alphabet=None, **kwargs)
