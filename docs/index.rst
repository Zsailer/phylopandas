.. phylopandas documentation master file, created by
   sphinx-quickstart on Mon Oct 30 16:22:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PhyloPandas
===========

**Bringing the Pandas DataFrame to phylogenetics.**

PhyloPandas provides a Pandas-like interface for reading various sequence formats into DataFrames. This enables easy manipulation of phylogenetic data using familiar Python/Pandas functions. Finally, phylogenetics for humans!

.. image:: _images/jlab.png
  :align: center

How does it work?
-----------------

Don't worry, we didn't reinvent the wheel. **PhyloPandas** is simply a DataFrame_
(great for human-accessible data storage) interface on top of Biopython_ (great for parsing/writing sequence data).

.. _DataFrame: https://github.com/pandas-dev/pandas
.. _Biopython: https://github.com/biopython/biopython

Basic Usage
-----------

1. Read any format:

.. code-block:: python

  import phylopandas as pd

  df1 = pd.read_fasta('sequences.fasta')
  df2 = pd.read_phylip('sequences.phy')

2. Write any format:

.. code-block:: python

  df1.to_clustal('sequences.clustal')

3. Convert formats:

.. code-block:: python


  df = phypd.read_fasta('sequences.fasta')
  df.to_phylip('sequences.phy')


4. Merge two **ordered** sequence files (like raw sequence file and its alignment).

.. code-block:: python

  # Read sequence file into dataframe
  df = pd.read_fasta('sequences.fasta')

  # Read alignment into dataframe
  align = pd.read_fasta('alignment.fasta')

  # Add alignment using standard pandas functions
  # NOTE: this assumes the alignment and sequence
  #       file are ordered.
  df = df.assign(alignment=align['sequence'])

5. Write out alignment in last example.

.. code-block:: python

  df.to_fasta('new_alignment.fasta', sequence_col='alignment')

Contributing
------------

It's *easy* to create new read/write functions and methods for PhyloPandas. If you
have a format you'd like to add, please submit PRs! There are many more formats
in Biopython that I haven't had the time to add myself, so please don't be afraid
to add then yourself! I thank you ahead of time!

Dependencies
------------

* BioPython
* Pandas

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   _pages/cookbook


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
