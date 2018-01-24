Cookbook
========

Merge two sequence files
------------------------

Use ``pandas.concat`` to merge two sequence files.

.. code-block:: python

  import phylopandas as ph

  seq1 = ph.read_fasta('seq1.fasta')
  seq2 = ph.read_fasta('seq2.fasta')

  # Merge two files
  seqs = seq1.concat(seq2, ignore_index=False)

  # Write to file.
  seqs.to_fasta('seqs.fasta')



Merge alignment and sequence data
---------------------------------

Add alignment column to sequence DataFrame.

.. code-block:: python

  import phylopandas as ph

  # Read sequences and alignments.
  seq = ph.read_fasta('sequences.fasta')
  ali = ph.read_fasta('alignment.fasta')

  # Merge data.
  seq.merge(ali, on='id')
