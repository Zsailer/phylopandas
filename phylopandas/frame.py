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


class DataFrame(pd.DataFrame):
    """ Two-dimensional size-mutable, potentially heterogeneous tabular data
    structure with labeled axes (rows and columns). Arithmetic operations
    align on both row and column labels. Can be thought of as a dict-like
    container for Series objects. The primary pandas data structure
    Parameters
    ----------
    data : numpy ndarray (structured or homogeneous), dict, or DataFrame
        Dict can contain Series, arrays, constants, or list-like objects
    index : Index or array-like
        Index to use for resulting frame. Will default to np.arange(n) if
        no indexing information part of input data and no index provided
    columns : Index or array-like
        Column labels to use for resulting frame. Will default to
        np.arange(n) if no column labels are provided
    dtype : dtype, default None
        Data type to force. Only a single dtype is allowed. If None, infer
    copy : boolean, default False
        Copy data from inputs. Only affects DataFrame / 2d ndarray input

    """
    def __init__(self, data=None, *args, **kwargs):
        # Prevent type errors thrown by pandas.DataFrame constructor.
        # Get core data.
        if isinstance(data, DataFrame) or isinstance(data, pd.DataFrame):
            data = data._data

        # Initialize with pandas init function.
        super(DataFrame, self).__init__(data=data, *args, **kwargs)

    def to_fasta(self, filename=None, sequence_col='sequence',
                 id_col='id', id_only=False):
        """Write to fasta format."""
        return _write(self, filename=filename, schema="fasta",
                      id_col=id_col, id_only=id_only,
                      sequence_col=sequence_col)

    def to_phylip(self, filename=None, sequence_col='sequence', id_col='id'):
        """Write to phylip format."""
        return _write(self, filename=filename, schema="phylip",
                      id_col=id_col, sequence_col=sequence_col)

    def to_clustal(self, filename=None, sequence_col='sequence'):
        """Write to alignment format of Clustal X and Clustal W."""
        return _write(self, filename=filename, schema="clustal",
                      sequence_col=sequence_col)

    def to_embl(self, alphabet, filename=None, sequence_col='sequence'):
        """Write to the EMBL flat file format."""
        return _write(self, alphabet=alphabet, filename=filename,
                      schema="embl", sequence_col=sequence_col)

    def to_nexus(self, alphabet, filename=None, sequence_col='sequence'):
        """Write to the NEXUS multiple alignment format."""
        return _write(self, alphabet=alphabet, filename=filename,
                      schema="nexus", sequence_col=sequence_col)

    def to_swiss(self, filename=None, sequence_col='sequence'):
        """Write Swiss-Prot aka UniProt format."""
        return _write(self, filename=filename, schema="swiss",
                      sequence_col=sequence_col)
