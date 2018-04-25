__doc__ = """
Functions for write sequence data to sequence files.
"""
import pandas as pd

# Import Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Alphabet


def _write_doc_template(schema):
    s = """Write to {} format.

    Parameters
    ----------
    filename : str
        File to write {} string to. If no filename is given, a {} string
        will be returned.

    sequence_col : str (default='sequence')
        Sequence column name in DataFrame.

    id_col : str (default='id')
        ID column name in DataFrame

    id_only : bool (default=False)
        If True, use only the ID column to label sequences in fasta.
    """.format(schema, schema, schema)
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

def _write_method(schema):
    """Add a write method for named schema to a class.
    """
    def method(
        self,
        filename=None,
        sequence_col='sequence',
        id_col='id',
        id_only=False,
        alphabet=None,
        **kwargs):
        # Use generic write class to write data.
        return _write(
            self._data,
            filename=filename,
            schema=schema,
            sequence_col=sequence_col,
            id_col=id_col,
            id_only=id_only,
            alphabet=alphabet,
            **kwargs
        )
    # Update docs
    method.__doc__ = _write_doc_template(schema)
    return method


def _write_function(schema):
    """Add a write method for named schema to a class.
    """
    def func(
        data,
        filename=None,
        sequence_col='sequence',
        id_col='id',
        id_only=False,
        alphabet=None,
        **kwargs):
        # Use generic write class to write data.
        return _write(
            data,
            filename=filename,
            schema=schema,
            sequence_col=sequence_col,
            id_col=id_col,
            id_only=id_only,
            alphabet=alphabet,
            **kwargs
        )
    # Update docs
    func.__doc__ = _write_doc_template(schema)
    return func


# Write functions to various formats.
to_fasta = _write_function('fasta')
to_phylip = _write_function('phylip')
to_clustal = _write_function('clustal')
to_embl = _write_function('embl')
to_nexus = _write_function('nexus')
to_swiss = _write_function('swiss')
to_fastq = _write_function('fastq')
