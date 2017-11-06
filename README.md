<img src="docs/_logo/banner.png">

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/phylopandas/Lobby)
[![Documentation Status](http://readthedocs.org/projects/phylopandas/badge/?version=latest)](http://phylopandas.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/Zsailer/phylopandas.svg?branch=master)](https://travis-ci.org/Zsailer/phylopandas)

**Bringing the [Pandas](https://github.com/pandas-dev/pandas) `DataFrame` to phylogenetics.**


PhyloPandas provides a Pandas-like interface for reading various sequence formats into DataFrames. This enables easy manipulation of phylogenetic data using familiar Python/Pandas functions. Finally, phylogenetics for humans!

<img src='docs/_images/jlab.png' align="middle">

## How does it work?

Don't worry, we didn't reinvent the wheel. **PhyloPandas** is simply a [DataFrame](https://github.com/pandas-dev/pandas) 
(great for human-accessible data storage) interface on top of [Biopython](https://github.com/biopython/biopython) (great for parsing/writing sequence data). 

When you import PhyloPandas, you import Pandas with a PhyloPandas flavor. That means, the usual `read_` functions 
are available ('read_csv', 'read_excel', etc.), but the returned DataFrame includes extra `to_` methods (`to_fasta`, `to_phylip`, etc.) 

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

If you have ideas for the project, please share them on the project's [Gitter chat](https://gitter.im/phylopandas/Lobby). 

It's *easy* to create new read/write functions and methods for PhyloPandas. If you 
have a format you'd like to add, please submit PRs! There are many more formats 
in Biopython that I haven't had the time to add myself, so please don't be afraid
to add them! I thank you ahead of time!

## Testing

PhyloPandas includes a small [pytest](https://docs.pytest.org/en/latest/) suite. Run these tests from base directory.
```
$ cd phylopandas
$ pytest
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

* [BioPython](https://github.com/biopython/biopython): Library for managing and manipulating biological data.
* [Pandas](https://github.com/pandas-dev/pandas): Flexible and powerful data analysis / manipulation library for Python
