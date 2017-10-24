# PhyloPandas: Pandas for phylogenetics

PhyloPandas brings the Pandas `DataFrame` to phylogenetics. 

It provides `read_` methods from a large list of sequence formats and
a subclass of Panda's `DataFrame` for easy manipulation of phylogenetic data.
Boom! Phylogenetics for humans. 

How does it work? There is no reinventing the wheel here. PhyloPandas simply provides a
bridge between BioPython (great for parsing/writing sequence data) and Pandas 
(great for human-accessible data storage).   

# Basic Usage

Read from any format:

```python
import phylopandas as phypd

df = phypd.read_fasta('sequences.fasta')
df = phypd.read_phylip('sequences.phy')

```

Convert formats:

```python
df = phypd.read_fasta('sequences.fasta')
df.to_phylip('sequences.phy')
```

# Install

Install from source:
```
git clone https://github.com/Zsailer/phylopandas
cd phylopandas
pip install -e .
```

# Dependencies

* BioPython
* Pandas
