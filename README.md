# PhyloPandas # 

**Bringing the [Pandas](https://github.com/pandas-dev/pandas) `DataFrame` to phylogenetics.**

PhyloPandas provides a Pandas-like interface for reading various sequence formats into DataFrames. This enables easy manipulation of phylogenetic data using familiar Python/Pandas functions. Finally, phylogenetics for humans!

## How does it work?

Don't worry, we didn't reinvent the wheel. **PhyloPandas** is simply places a [DataFrame](https://github.com/pandas-dev/pandas) 
(great for human-accessible data storage) interface on top of [Biopython's SeqIO](https://github.com/biopython/biopython) module  (great for parsing/writing sequence data). 


## Basic Usage

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
