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


def pandas_df_to_biopython_seqrecord(
    df,
    id_col='uid',
    sequence_col='sequence',
    extra_data=None,
    alphabet=None,
    ):
    """Convert pandas dataframe to biopython seqrecord for easy writing.

    Parameters
    ----------
    df : Dataframe
        Pandas dataframe to convert

    id_col : str
        column in dataframe to use as sequence label

    sequence_col str:
        column in dataframe to use as sequence data

    extra_data : list
        extra columns to use in sequence description line

    alphabet :
        biopython Alphabet object

    Returns
    -------
    seq_records :
        List of biopython seqrecords.
    """
    seq_records = []

    for i, row in df.iterrows():
        # Tries getting sequence data. If a TypeError at the seqrecord
        # creation is thrown, it is assumed that this row does not contain
        # sequence data and therefore the row is ignored.
        try:
            # Get sequence
            seq = Seq(row[sequence_col], alphabet=alphabet)

            # Get id
            id = row[id_col]

            # Build a description
            description = ""
            if extra_data is not None:
                description = " ".join([row[key] for key in extra_data])

            # Build a record
            record = SeqRecord(
                seq=seq,
                id=id,
                description=description,
            )
            seq_records.append(record)
        except TypeError:
            pass

    return seq_records

def pandas_series_to_biopython_seqrecord(
    series,
    id_col='uid',
    sequence_col='sequence',
    extra_data=None,
    alphabet=None
    ):
    """Convert pandas series to biopython seqrecord for easy writing.

    Parameters
    ----------
    series : Series
        Pandas series to convert

    id_col : str
        column in dataframe to use as sequence label

    sequence_col : str
        column in dataframe to use as sequence data

    extra_data : list
        extra columns to use in sequence description line

    Returns
    -------
    seq_records :
        List of biopython seqrecords.
    """
    # Get sequence
    seq = Seq(series[sequence_col], alphabet=alphabet)

    # Get id
    id = series[id_col]

    # Build a description
    description = ""
    if extra_data is not None:
        description = " ".join([series[key] for key in extra_data])

    # Build a record
    record = SeqRecord(
        seq=seq,
        id=id,
        description=description,
    )

    seq_records = [record]
    return seq_records

def _write(
    data,
    filename=None,
    schema='fasta',
    id_col='uid',
    sequence_col='sequence',
    extra_data=None,
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

    # Build a list of records from a pandas DataFrame
    if type(data) is pd.DataFrame:
        seq_records = pandas_df_to_biopython_seqrecord(
            data,
            id_col=id_col,
            sequence_col=sequence_col,
            extra_data=extra_data,
            alphabet=alphabet,
        )

    # Build a record from a pandas Series
    elif type(data) is pd.Series:
        seq_records = pandas_series_to_biopython_seqrecord(
            data,
            id_col=id_col,
            sequence_col=sequence_col,
            extra_data=extra_data,
            alphabet=alphabet,
        )

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
        schema=schema,
        id_col='uid',
        sequence_col='sequence',
        extra_data=None,
        alphabet=None,
        **kwargs):
        # Use generic write class to write data.
        return _write(
            self._data,
            filename=filename,
            schema=schema,
            id_col=id_col,
            sequence_col=sequence_col,
            extra_data=extra_data,
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
        schema=schema,
        id_col='uid',
        sequence_col='sequence',
        extra_data=None,
        alphabet=None,
        **kwargs):
        # Use generic write class to write data.
        return _write(
            data,
            filename=filename,
            schema=schema,
            id_col=id_col,
            sequence_col=sequence_col,
            extra_data=extra_data,
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
to_nexus_seq = _write_function('nexus')
to_swiss = _write_function('swiss')
to_fastq = _write_function('fastq')
