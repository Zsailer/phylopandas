The PhyloPandas DataFrame
=========================

The phylopandas dataframe is the core datastructure in this package. It defines
a set of columns (or grammar) for phylogenetic data. A few advantages of
defining such a grammar is: 1) we can leverage powerful+interative
visualization tools like Vega and 2) we standardize phylogenetic data in a
familiar format.

Columns of a Phylopandas DataFrame
----------------------------------

When reading sequence data, the following information will be stored on the dataframe.

1. ``sequence`` : DNA or protein sequence.
2. ``id``: user defined label or identifier.
3. ``description``: user defined description.

When reading tree data, the following information will be stored on the dataframe.

1. ``type`` : label describing the type of node; either "leaf" or "node".
2. ``parent`` : label of parent node.
3. ``branch_length`` : distance from parent node.

PhyloPandas indexes each sequence using a randomly generated 10 character key.

If reading tree data from a PhyloPandas DataFrame containing sequence data, the
two dataframes will be merged on the randomly generated index (unless otherwise specified).

If reading sequence data from a PhyloPandas DataFrae containing tree data, the two dataframes will be merged on the randomly generated index (unless otherwise specified).
