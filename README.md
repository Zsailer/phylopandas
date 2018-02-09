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

PhyloPandas does two things:
1. It offers new `read` functions to read sequence data directory into a DataFrame.
2. It attaches a new `phylo` **accessor** to the Pandas DataFrame. This accessor provides writing methods for sequencing data (powered by Biopython).

## Basic Usage

Read in a sequence file.
```python
import phylopandas as ph

df1 = ph.read_fasta('sequences.fasta')
df2 = ph.read_phylip('sequences.phy')
```

Write to various sequence file formats.

```python
df1.phylo.to_clustal('sequences.clustal')
```

Convert between formats.

```python
# Read a format.
df = ph.read_fasta('sequences.fasta')

# Write to a different format.
df.phylo.to_phylip('sequences.phy')
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

- [BioPython](https://github.com/biopython/biopython): Library for managing and manipulating biological data.
- [Pandas](https://github.com/pandas-dev/pandas): Flexible and powerful data analysis / manipulation library for Python
- [pandas_flavor](https://github.com/Zsailer/pandas_flavor): Flavor pandas objects with new accessors using pandas' new register API (with backwards compatibility).
