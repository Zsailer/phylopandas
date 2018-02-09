__doc__ = """
Functions for write sequence data to sequence files.
"""
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

def _write(dataframe, filename=None, schema='fasta', sequence_col='sequence',
           id_col='id', id_only=False, alphabet=None, **kwargs):
    """Write a PhyloPandas DataFrame to a sequence file.

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
    # Check Alphabet if given
    if alphabet is None:
        alphabet = Bio.Alphabet.Alphabet()

    elif alphabet in ['dna', 'rna', 'protein', 'nucleotide']:
        alphabet = getattr(Bio.Alphabet, 'generic_{}'.format(alphabet))

    else:
        raise Exception(
            "The alphabet is not recognized. Must be 'dna', 'rna', "
            "'nucleotide', or 'protein'.")

    # Iterate through dataframe rows
    seq_records = []
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
