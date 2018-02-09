# Import pandas
from pandas_flavor import register_dataframe_accessor

from functools import wraps
from . import seqio
from . import treeio


@register_dataframe_accessor('phylo')
class PhyloPandasMethods(object):
    """PhyloPandas accessor to the Pandas DataFrame.

    This accessor adds reading/writing methods to the pandas DataFrame that
    are specific to phylogenetic data.
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

    @wraps(seqio.write.to_nexus)
    def to_nexus(self, *args, **kwargs):
        return seqio.write.to_nexus(self._data, *args, **kwargs)

    @wraps(seqio.write.to_swiss)
    def to_swiss(self, *args, **kwargs):
        return seqio.write.to_swiss(self._data, *args, **kwargs)

    @wraps(seqio.write.to_nexus)
    def to_nexus(self, *args, **kwargs):
        return seqio.write.to_nexus(self._data, *args, **kwargs)

    @wraps(seqio.write.to_fastq)
    def to_fastq(self, *args, **kwargs):
        return seqio.write.to_fastq(self._data, *args, **kwargs)
