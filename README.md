# PhyloPandas # 

**Bringing the [Pandas](https://github.com/pandas-dev/pandas) `DataFrame` to phylogenetics.**

PhyloPandas provides a Pandas-like interface for reading various sequence formats into DataFrames. This enables easy manipulation of phylogenetic data using familiar Python/Pandas functions. Finally, phylogenetics for humans!

<img src='docs/_images/jlab.png' align="middle">

## How does it work?

Don't worry, we didn't reinvent the wheel. **PhyloPandas** is simply a [DataFrame](https://github.com/pandas-dev/pandas) 
(great for human-accessible data storage) interface on top of [Biopython](https://github.com/biopython/biopython) (great for parsing/writing sequence data). 

## Basic Usage

1. Read any format:
```python
import phylopandas as pd

df1 = pd.read_fasta('sequences.fasta')
df2 = pd.read_phylip('sequences.phy')
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
df = pd.read_fasta('sequences.fasta')

# Read alignment into dataframe
align = pd.read_fasta('alignment.fasta')

# Add alignment using standard pandas functions
# NOTE: this assumes the alignment and sequence
#       file are ordered.
df = df.assign(alignment=align['sequence'])
```
5. Write out alignment in last example.
```python
df.to_fasta('new_alignment.fasta', sequence_col='alignment')
``` 

## Contributing

It's *easy* to create new read/write functions and methods for PhyloPandas. If you 
have a format you'd like to add, please submit PRs! There are many more formats 
in Biopython that I haven't had the time to add myself, so please don't be afraid
to add then yourself! I thank you ahead of time!

## Install

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

## Dependencies

* BioPython
* Pandas
