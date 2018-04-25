# Import pandas
from pandas_flavor import register_dataframe_accessor, register_series_accessor

from functools import wraps
from . import seqio
from . import treeio


def verify_phylopandas_function(f):
    """"""
    @wraps(f)
    def inner(data, *args, **kwargs):
        # Sanity check.
        if not hasattr(data, 'phylo'):
            raise Exception("Object is not a PhyloPandas dataframe.")
        return f(args, kwargs)
    return inner


@register_series_accessor('phylo')
class PhyloPandasSeriesMethods(object):
    """
    """
    def __init__(self, data):
        self._data = data

    @wraps(seqio.write.to_fasta)
    def to_fasta(self, *args, **kwargs):
        return seqio.write.to_fasta(self._data, *args, **kwargs)

    @wraps(seqio.write.to_phylip)
    def to_phylip(self, *args, **kwargs):
        return seqio.write.to_phylip(self._data, *args, **kwargs)

    @wraps(seqio.write.to_clustal)
    def to_clustal(self, *args, **kwargs):
        return seqio.write.to_clustal(self._data, *args, **kwargs)

    @wraps(seqio.write.to_embl)
    def to_embl(self, *args, **kwargs):
        return seqio.write.to_embl(self._data, *args, **kwargs)

    @wraps(seqio.write.to_swiss)
    def to_swiss(self, *args, **kwargs):
        return seqio.write.to_swiss(self._data, *args, **kwargs)

    @wraps(seqio.write.to_nexus)
    def to_nexus(self, *args, **kwargs):
        return seqio.write.to_nexus(self._data, *args, **kwargs)

    @wraps(seqio.write.to_fastq)
    def to_fastq(self, *args, **kwargs):
        return seqio.write.to_fastq(self._data, *args, **kwargs)


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

    @wraps(seqio.write.to_fasta)
    def to_fasta(self, *args, **kwargs):
        return seqio.write.to_fasta(self._data, *args, **kwargs)

    @wraps(seqio.write.to_phylip)
    def to_phylip(self, *args, **kwargs):
        return seqio.write.to_phylip(self._data, *args, **kwargs)

    @wraps(seqio.write.to_clustal)
    def to_clustal(self, *args, **kwargs):
        return seqio.write.to_clustal(self._data, *args, **kwargs)

    @wraps(seqio.write.to_embl)
    def to_embl(self, *args, **kwargs):
        return seqio.write.to_embl(self._data, *args, **kwargs)

    @wraps(seqio.write.to_swiss)
    def to_swiss(self, *args, **kwargs):
        return seqio.write.to_swiss(self._data, *args, **kwargs)

    @wraps(seqio.write.to_nexus)
    def to_nexus(self, *args, **kwargs):
        return seqio.write.to_nexus(self._data, *args, **kwargs)

    @wraps(seqio.write.to_fastq)
    def to_fastq(self, *args, **kwargs):
        return seqio.write.to_fastq(self._data, *args, **kwargs)

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
