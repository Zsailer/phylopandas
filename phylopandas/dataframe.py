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

def _write(dataframe, filename, format, sequence_col='sequence', **kwargs):
    """Write a PhyloPandas DataFrame to a sequence file.
    """
    seq_records = []
    # Iterate through dataframe rows
    for i, row in dataframe.iterrows():
        
        # Create a SeqRecord and append to list.
        seq_records.append(SeqRecord(
            Seq(row[sequence_col]),
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

def read_paml(filename, **kwargs):
    """Read PAML output and get ancestors as DataFrame. Returns DataFrame and tree as a string.
    """
    # Read paml output.
    with open(filename, 'r') as f:
        data = f.read()
    
    # Match tree with ancestor labels
    regex = re.compile("tree with node labels for Rod Page's TreeView\n.+\n")
    match = regex.search(data)
    tree = match.group().split("\n")[1]

    # Initialize the ancestor_data dictionary
    ancestor_data = {}

    # Compile a regular expression to find blocks of data for internal nodes
    node_regex = re.compile("""Prob distribution at node [0-9]+, by site[-\w():.\s]+\n\n""")
    # Strip the node number from this block of data.
    node_num_regex = re.compile("[0-9]+")

    # Iterate through each block of internal node data
    index, sequences, posteriors = [], [], []
    for node in node_regex.findall(data):
        # Initialize a dictionary for site data
        site_data = {}

        # Compile regex for matching site data
        site_regex = re.compile("(?:\w\(\w.\w{3}\) )+")

        site_num = 0
        # Iterate through each match for site data.
        seq, post = [], []
        for site in site_regex.findall(node):
            # Iterate through residues
            scores = [float(site[i+2:i+7]) for i in range(0,len(site), 9)]
            j, p = max(enumerate(scores), key=lambda item: item[1])
            seq.append(site[j*9])
            post.append(p)

        # Add site data to ancestor_data
        index.append(int(node_num_regex.search(node).group(0))) # Get node_number defined by PAML
        sequences.append("".join(seq))
        posteriors.append(sum(post)/len(post))
        
    df = DataFrame({"sequence":sequences, "posterior":posteriors}, index=index)
    return tree, df


class DataFrame(pd.DataFrame):
    
    def to_fasta(self, f, sequence_col='sequence'):
        """Write to fasta format."""
        _write(self, f, format="fasta")
        
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
        
    
    
