# PhyloPandas: Pandas for phylogenetics

Phylogenetics is cursed by a myriad of file formats. It seems like each piece of 
phylogenetic software defines its own format and won't accept extra metadata maybe
useful for interpreting output. I found myself constantly writing parsers and writers, 
stripping out metadata, shortening labels, and injecting metadata at the end.
To make things worse, most of these formats are unreadable to any human. 

PhyloPandas simplifies this problem. It uses `pandas.DataFrame` to centralize
the data. PhyloPandas subclasses `DataFrame` and appends a bunch of methods for
writing to phylogenetic formats necessary for working with phylogenetic software.
Then, it seamless reads the output and adds new columns of data to the DataFrame.

When you're finished, write your DataFrame to a CSV or Excel file using the usual
pandas interface. Boom! Human readable phylogenetics. 

How does it work? There is no reinventing the wheel here. PhyloPandas simply bridges
BioPython (great for parsing sequence data) and Pandas.  

# Basic Usage

```python
import phylopandas as phypd

df = phypd.read_fasta('sequences.fasta')
df
```

Write the dataframe to Phylip format.

```python
df.to_phylip('sequences.phy')
```
