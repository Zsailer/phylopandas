from Bio import SeqIO

from .dataframe import DataFrame

def _read(filename, format, **kwargs):
    """Use BioPython's sequence parsing module to convert any file format to 
    a Pandas DataFrame.
    """
    # Prepare DataFrame fields.
    data = {'id':[], 'sequence':[], 'description':[], 'name':[], 'db_id':[]}
    
    # Parse Fasta file.
    for i, s in enumerate(SeqIO.parse(filename, 'fasta', **kwargs)):
        data['db_id'].append(s.id)
        data['sequence'].append(str(s.seq))
        data['description'].append(s.description)
        data['name'].append(s.name)
        data['id'].append('XX{:08d}'.format(i))
    
    # Port to DataFrame.
    return DataFrame(data)

def read_fasta(filename, **kwargs):
    """
    """
    return _read(filename, format='fasta', **kwargs)

def read_phylip():
    """
    """
    return _read(filename, format='fasta', **kwargs)
