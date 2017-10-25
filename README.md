# PhyloPandas

*Bringing [Pandas](https://github.com/pandas-dev/pandas) `DataFrame` to phylogenetics.*

Read sequence formats into Pandas `DataFrame` for easy manipulation of phylogenetic data. **Finally, phylogenetics for humans!**

# how does it work?

Don't worry, we didn't reinvent the wheel here. PhyloPandas simply bridges [BioPython](https://github.com/biopython/biopython) (great for parsing/writing sequence data) and [Pandas](https://github.com/pandas-dev/pandas) 
(great for human-accessible data storage).   

# things you can do

1. Read any format:
```python
import phylopandas as phypd

df1 = phypd.read_fasta('sequences.fasta')
df2 = phypd.read_phylip('sequences.phy')
```
2. Write any format:
```python
df1.to_clustal('sequences.clustal')
```
3. Convert formats:
```python
df = phypd.read_fasta('sequences.fasta')
df.to_phylip('sequences.phy')
```
4. Merge two **ordered** sequence files (like raw sequence file and its alignment).
```python
# Read sequence file into dataframe
df = phypd.read_fasta('sequences.fasta')

# Read alignment into dataframe
align = phypd.read_fasta('alignment.fasta')

# Add alignment using standard pandas functions
# NOTE: this assumes the alignment and sequence
#       file are ordered.
df = df.assign(alignment=align['sequence'])
```
5. Write out alignment in last example.
```python
df.to_fasta('new_alignment.fasta', sequence_col='alignment')
``` 

# installation

Install from PyPi:

```
pip install phylopandas
```

Install from source:

```
git clone https://github.com/Zsailer/phylopandas
cd phylopandas
pip install -e .
```

# dependencies

* BioPython
* Pandas
