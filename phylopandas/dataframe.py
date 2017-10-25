import os
import re 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd

def _read(filename, format, seq_label='sequence', **kwargs):
    """Use BioPython's sequence parsing module to convert any file format to 
    a Pandas DataFrame.
    """
    # Prepare DataFrame fields.
    data = {'id':[], seq_label:[], 'description':[], 'name':[], 'db_id':[]}
    
    # Parse Fasta file.
    for i, s in enumerate(SeqIO.parse(filename, format=format, **kwargs)):
        data['db_id'].append(s.id)
        data[seq_label].append(str(s.seq))
        data['description'].append(s.description)
        data['name'].append(s.name)
        data['id'].append('XX{:08d}'.format(i))
    
    # Port to DataFrame.
    return DataFrame(data)

def _write(dataframe, filename, format, sequence_col='sequence', id_only=False, **kwargs):
    """Write a PhyloPandas DataFrame to a sequence file.
    """
    seq_records = []
    # Iterate through dataframe rows
    for i, row in dataframe.iterrows():
        seq = Seq(row[sequence_col])
        # Create a SeqRecord and append to list.
        if id_only:
            seq_records.append(SeqRecord(
                seq,
                id=row['id'],
                name='',
                description=''))
        else:
            seq_records.append(SeqRecord(
                seq,
                id=row['id'],
                name=row['name'],
                description=row['description']))
            
    # Write out sequences        
    SeqIO.write(seq_records, filename, format=format, **kwargs)

def read_fasta(filename, **kwargs):
    """Read fasta format."""
    return _read(filename, format='fasta', **kwargs)

def read_phylip(filename, **kwargs):
    """Read phylip format."""
    return _read(filename, format='phylip', **kwargs)

def read_clustal(filename, **kwargs):
    """Read clustal format."""
    return _read(filename, format='clustal', **kwargs)

def read_embl(filename, **kwargs):
    """Read the EMBL flat file format."""
    return _read(filename, format='embl', **kwargs)

def read_nexus(filename, **kwargs):
    """Read the EMBL flat file format."""
    return _read(filename, format='nexus', **kwargs)

def read_swiss(filename, **kwargs):
    """Read Swiss-Prot aka UniProt format."""
    return _read(filename, format='nexus', **kwargs)


class DataFrame(pd.DataFrame):
    
    def to_fasta(self, f, sequence_col='sequence', id_only=True):
        """Write to fasta format."""
        _write(self, f, format="fasta", id_only=True)
        
    def to_phylip(self, f, sequence_col='sequence'):
        """Write to phylip format."""
        _write(self, f, format="phylip")

    def to_clustal(self, f, sequence_col='sequence'):
        """Write to alignment format of Clustal X and Clustal W."""
        _write(self, f, format="clustal")
        
    def to_embl(self, f, sequence_col='sequence'):
        """Write to the EMBL flat file format."""
        _write(self, f, format="embl")

    def to_nexus(self, f, sequence_col='sequence'):
        """Write to the NEXUS multiple alignment format."""
        _write(self, f, format="nexus")
        
    def to_swiss(self, f, sequence_col='sequence'):
        """Write Swiss-Prot aka UniProt format."""
        _write(self, f, format="swiss")
        
    
    
