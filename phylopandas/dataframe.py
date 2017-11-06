import os
import re 
from functools import wraps

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import Bio.Alphabet

import pandas as pd

def _read(filename, schema, seq_label='sequence', alphabet=None, **kwargs):
    """Use BioPython's sequence parsing module to convert any file format to 
    a Pandas DataFrame.
    """
    # Check Alphabet if given
    if alphabet is None: 
        alphabet = Bio.Alphabet.Alphabet()
    elif alphabet in ['dna', 'rna', 'protein', 'nucleotide']:
        alphabet = getattr(Bio.Alphabet, 'generic_{}'.format(alphabet))
    else:
        raise Exception("The alphabet is not recognized. Must be 'dna', 'rna', 'nucleotide', or 'protein'.")
    
    kwargs.update(alphabet=alphabet)

    # Prepare DataFrame fields.
    data = {'id':[], seq_label:[], 'description':[], 'name':[]}
    
    # Parse Fasta file.
    for i, s in enumerate(SeqIO.parse(filename, format=schema, **kwargs)):
        data['id'].append(s.id)
        data[seq_label].append(str(s.seq))
        data['description'].append(s.description)
        data['name'].append(s.name)
    
    # Port to DataFrame.
    return DataFrame(data)

def _write(dataframe, filename=None, schema='fasta', sequence_col='sequence', id_col='id', id_only=False, alphabet=None, **kwargs):
    """Write a PhyloPandas DataFrame to a sequence file.
    """
    # Check Alphabet if given
    if alphabet is None: 
        alphabet = Bio.Alphabet.Alphabet()
    elif alphabet in ['dna', 'rna', 'protein', 'nucleotide']:
        alphabet = getattr(Bio.Alphabet, 'generic_{}'.format(alphabet))
    else:
        raise Exception("The alphabet is not recognized. Must be 'dna', 'rna', 'nucleotide', or 'protein'.")
    
    seq_records = []
    # Iterate through dataframe rows
    for i, row in dataframe.iterrows():
        seq = Seq(row[sequence_col], alphabet)
        # Create a SeqRecord and append to list.
        if id_only:
            record = SeqRecord(seq, id=row[id_col], name='', description='')
        else:
            record = SeqRecord(seq, id=row[id_col], name=row['name'], description=row['description'])
        seq_records.append(record)
            
    # Write to disk or return string
    if filename != None:
        SeqIO.write(seq_records, filename, format=schema, **kwargs)
    else:
        return "".join([s.format(schema) for s in seq_records])

def read_fasta(filename, **kwargs):
    """Read fasta format."""
    return _read(filename, schema='fasta', **kwargs)

def read_phylip(filename, **kwargs):
    """Read phylip format."""
    return _read(filename, schema='phylip', **kwargs)

def read_clustal(filename, **kwargs):
    """Read clustal format."""
    return _read(filename, schema='clustal', **kwargs)

def read_embl(filename, **kwargs):
    """Read the EMBL flat file format."""
    return _read(filename, schema='embl', **kwargs)

def read_nexus(filename, **kwargs):
    """Read the EMBL flat file format."""
    return _read(filename, schema='nexus', **kwargs)

def read_swiss(filename, **kwargs):
    """Read Swiss-Prot aka UniProt format."""
    return _read(filename, schema='nexus', **kwargs)

def read_blast_xml(filename, **kwargs):
    """Read BLAST XML format."""
    # Read file.
    with open(filename, 'r') as f:
        blast_record = NCBIXML.read(f)    

    # Prepare DataFrame fields.
    data = {'accession':[], 
        'hit_def':[], 
        'hit_id':[], 
        'title':[],
        'length':[],
        'e_value':[],
        'sequence':[]}
    
    # Get alignments from blast result.
    for i, s in enumerate(blast_record.alignments):
        data['accession'] = s.accession
        data['hit_def'] = s.hit_def
        data['hit_id'] = s.hit_id
        data['title'] = s.title
        data['length'] = s.length
        data['e_value'] = s.hsps[0].expect
        data['sequence'] = s.hsps[0].sbjct
        
    # Port to DataFrame.
    return DataFrame(data)    


class DataFrame(pd.DataFrame):
    __doc__ == pd.DataFrame.__doc__
    
    def to_fasta(self, filename=None, sequence_col='sequence', id_col='id', id_only=False):
        """Write to fasta format."""
        return _write(self, filename=filename, schema="fasta", id_col=id_col, id_only=id_only, sequence_col=sequence_col)
        
    def to_phylip(self, filename=None, sequence_col='sequence', id_col='id'):
        """Write to phylip format."""
        return _write(self, filename=filename, schema="phylip", id_col=id_col, sequence_col=sequence_col)

    def to_clustal(self, filename=None, sequence_col='sequence'):
        """Write to alignment format of Clustal X and Clustal W."""
        return _write(self, filename=filename, schema="clustal", sequence_col=sequence_col)
        
    def to_embl(self, alphabet, filename=None, sequence_col='sequence'):
        """Write to the EMBL flat file format."""
        return _write(self, alphabet=alphabet, filename=filename, schema="embl", sequence_col=sequence_col)

    def to_nexus(self, alphabet, filename=None, sequence_col='sequence'):
        """Write to the NEXUS multiple alignment format."""
        return _write(self, alphabet=alphabet, filename=filename, schema="nexus", sequence_col=sequence_col)
        
    def to_swiss(self, filename=None, sequence_col='sequence'):
        """Write Swiss-Prot aka UniProt format."""
        return _write(self, filename=filename, schema="swiss", sequence_col=sequence_col)
