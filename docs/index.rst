.. phylopandas documentation master file, created by
   sphinx-quickstart on Mon Oct 30 16:22:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: _logo/banner.png


*Bringing the Pandas DataFrame to phylogenetics.*

PhyloPandas provides a Pandas-like interface for reading various sequence formats into DataFrames. This enables easy manipulation of phylogenetic data using familiar Python/Pandas functions. Finally, phylogenetics for humans!

.. image:: _images/jlab.png
  :align: center

|

How does it work?
-----------------

Don't worry, we didn't reinvent the wheel. **PhyloPandas** is simply a DataFrame_
(great for human-accessible data storage) interface on top of Biopython_ (great for parsing/writing sequence data).

.. _DataFrame: https://github.com/pandas-dev/pandas
.. _Biopython: https://github.com/biopython/biopython

Basic Usage
~~~~~~~~~~~

Read sequence file into DataFrame.

.. code-block:: python

  import phylopandas as ph

  df1 = ph.read_fasta('sequences.fasta')

Write ``phylopandas.DataFrame`` data to sequence file.

.. code-block:: python

  df1.to_clustal('sequences.clustal')

Convert between two sequence formats.

.. code-block:: python

  # Read from fasta.
  df = phypd.read_fasta('sequences.fasta')

  # Write to phylip.
  df.to_phylip('sequences.phy')

See the Cookbook_ page for more things you can do.

.. _Cookbook: _pages/cookbook.html

Contributing
~~~~~~~~~~~~

It's *easy* to create new read/write functions and methods for PhyloPandas. If you
have a format you'd like to add, please submit PRs! There are many more formats
in Biopython that I haven't had the time to add myself, so please don't be afraid
to add then yourself! I thank you ahead of time!


Table of Contents
~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   _pages/dataframe
   _pages/cookbook


Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
