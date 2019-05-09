# Import pandas
import pandas as pd
from pandas_flavor import register_dataframe_accessor, register_series_accessor

from . import seqio
from . import treeio


try: 
    from phylovega import TreeChart
except ImportError:
    TreeChart = None


@register_series_accessor('phylo')
class PhyloPandasSeriesMethods(object):
    """
    """
    def __init__(self, data):
        self._data = data

    # -----------------------------------------------------------
    # Extra write methods.
    # -----------------------------------------------------------

    to_fasta = seqio.write._write_method('fasta')
    to_phylip = seqio.write._write_method('phylip')
    to_clustal = seqio.write._write_method('clustal')
    to_embl = seqio.write._write_method('embl')
    to_nexus = seqio.write._write_method('nexus')
    to_swiss = seqio.write._write_method('swiss')
    to_fastq = seqio.write._write_method('fastq')
    to_fasta_twoline = seqio.write._write_method('fasta-2line')
    to_phylip_sequential = seqio.write._write_method('phylip-sequential')
    to_phylip_relaxed = seqio.write._write_method('phylip-relaxed')


@register_dataframe_accessor('phylo')
class PhyloPandasDataFrameMethods(object):
    """PhyloPandas accessor to the Pandas DataFrame.

    This accessor adds reading/writing methods to the pandas DataFrame that
    are specific to phylogenetic data.
    """
    def __init__(self, data):
        self._data = data

    # -----------------------------------------------------------
    # Extra read methods.
    # -----------------------------------------------------------

    # Sequence file reading methods
    read_fasta = seqio.read._read_method('fasta')
    read_phylip = seqio.read._read_method('phylip')
    read_clustal = seqio.read._read_method('clustal')
    read_embl = seqio.read._read_method('embl')
    read_nexus_seq = seqio.read._read_method('nexus')
    read_swiss = seqio.read._read_method('swiss')
    read_fastq = seqio.read._read_method('fastq')
    read_fasta_twoline = seqio.read._read_method('fasta-2line')
    read_phylip_sequential = seqio.read._read_method('phylip-sequential')
    read_phylip_relaxed = seqio.read._read_method('phylip-relaxed')

    # Tree file reading methods.
    read_newick = treeio.read._read_method('newick')
    read_nexus_tree = treeio.read._read_method('nexus')
    
    def read_dendropy(
        self,
        add_node_labels=True,
        combine_on='index',
        use_uids=True):
        df0 = self._data
        df1 = treeio.read.read_dendropy(
            self._data,
            add_node_labels=add_node_labels,
            use_uids=use_uids
        )
        return df0.phylo.combine(df1, on=combine_on)


    # -----------------------------------------------------------
    # Extra write methods.
    # -----------------------------------------------------------

    to_fasta = seqio.write._write_method('fasta')
    to_phylip = seqio.write._write_method('phylip')
    to_clustal = seqio.write._write_method('clustal')
    to_embl = seqio.write._write_method('embl')
    to_nexus_seq = seqio.write._write_method('nexus')
    to_swiss = seqio.write._write_method('swiss')
    to_fastq = seqio.write._write_method('fastq')
    to_fasta_twoline = seqio.write._write_method('fasta-2line')
    to_phylip_sequential = seqio.write._write_method('phylip-sequential')
    to_phylip_relaxed = seqio.write._write_method('phylip-relaxed')

    # Tree file reading methods.
    to_newick = treeio.write._write_method('newick')
    to_nexus_tree = treeio.write._write_method('nexus')

    def to_dendropy(
        self,
        taxon_col='uid',
        taxon_annotations=[],
        node_col='uid',
        node_annotations=[],
        branch_lengths=True):
        return treeio.write.to_dendropy(
            self._data,
            taxon_col=taxon_col,
            taxon_annotations=taxon_annotations,
            node_col=node_col,
            node_annotations=node_annotations,
            branch_lengths=branch_lengths,
        )

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


    def combine(self, other, on='index'):
        """Combine two dataframes. Update the first dataframe with second.
        New columns are added to the right of the first dataframe. Overlapping 
        columns update the values of the columns.

        Technical note: maintains order of columns, appending new dataframe to
        old.

        Parameters
        ----------
        other : DataFrame
            Index+Columns that match self will be updated with new values.
            New rows will be added separately.

        on : str
            Column to update index.
        """
        # Determine column labels for new dataframe (Maintain order of columns)
        column_idx = {k: None for k in self._data.columns}
        column_idx.update({k: None for k in other.columns})
        column_idx = list(column_idx.keys())

        df0 = self._data.copy()
        df1 = other.copy()

        # Set index to whatever column is given
        df0 = df0.set_index(on, inplace=False, drop=False)
        df1 = df1.set_index(on, inplace=False, drop=False)

        # Write out both dataframes to dictionaries
        data0 = df0.to_dict(orient="index")
        data1 = df1.to_dict(orient="index")

        # Update.
        for key in data1.keys():
            try:
                data0[key].update(data1[key])
            except KeyError:
                data0[key] = data1[key]

        # Build new dataframe
        df = pd.DataFrame(data0).T

        # Check for missing columns
        for key in column_idx:
            if key not in df.columns:
                df[key] = None

        # Reset the index.
        df.reset_index(inplace=True)

        # Return dataframe (maintaining original order)
        return df[column_idx]

    def show(self, **kwargs):
        __doc__ = TreeChart.__doc__
        # Show the tree using phylovega.
        try:
            return TreeChart(self._data.to_dict(orient='records'), **kwargs)
        except NameError:
            raise Exception("Looks like phylovega couldn't be imported. Is phylovega installed?")