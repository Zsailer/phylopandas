# Import pandas
from pandas_flavor import register_dataframe_accessor, register_series_accessor

from functools import wraps
from .seqio.write import _write_method
from . import seqio
from . import treeio


@register_series_accessor('phylo')
class PhyloPandasSeriesMethods(object):
    """
    """
    def __init__(self, data):
        self._data = data

    # -----------------------------------------------------------
    # Extra read/write methods.
    # -----------------------------------------------------------

    to_fasta = _write_method('fasta')
    to_phylip = _write_method('phylip')
    to_clustal = _write_method('clustal')
    to_embl = _write_method('embl')
    to_nexus = _write_method('nexus')
    to_swiss = _write_method('swiss')
    to_fastq = _write_method('fastq')


@register_dataframe_accessor('phylo')
class PhyloPandasDataFrameMethods(object):
    """PhyloPandas accessor to the Pandas DataFrame.

    This accessor adds reading/writing methods to the pandas DataFrame that
    are specific to phylogenetic data.
    """
    def __init__(self, data):
        self._data = data

    # -----------------------------------------------------------
    # Extra read/write methods.
    # -----------------------------------------------------------

    to_fasta = _write_method('fasta')
    to_phylip = _write_method('phylip')
    to_clustal = _write_method('clustal')
    to_embl = _write_method('embl')
    to_nexus = _write_method('nexus')
    to_swiss = _write_method('swiss')
    to_fastq = _write_method('fastq')

    # -----------------------------------------------------------
    # Useful dataframe methods specific to sequencing data.
    # -----------------------------------------------------------

    def match_value(self, column, value):
        """Return a subset dataframe that column values match the given value.

        Parameters
        ----------
        column : string
            column to search for matches

        value : float, int, list, etc.
            values to match.
        """
        # Get column
        col = self._data[column]

        # Get items in a list
        try:
            idx = col[col.isin(value)].index

        # Or value is a single item?
        except TypeError:
            idx = col[col == value].index

        return self._data.loc[idx]
