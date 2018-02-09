__doc__ = """
PhyloPandas
===========

*Pandas DataFrames for phylogenetics.*

PhyloPandas provides a Pandas-like interface for reading various sequence
formats into DataFrames. This enables easy manipulation of phylogenetic data
using familiar Python/Pandas functions. Finally, phylogenetics for humans!


How does it work?
-----------------

Don't worry, we didn't reinvent the wheel. PhyloPandas is simply a DataFrame
(great for human-accessible data storage) interface on top of Biopython
(great for parsing/writing sequence data).

When you import PhyloPandas, you import Pandas with a PhyloPandas flavor.
That means, the usual read_ functions are available ('read_csv',
'read_excel', etc.), but the returned DataFrame includes extra to_ methods
(to_fasta, to_phylip, etc.)
"""
# Import new read functions
from pandas import DataFrame

# Register PhyloPandas Methods
from .core import PhyloPandasMethods as _
from .seqio import *
from .treeio import *
