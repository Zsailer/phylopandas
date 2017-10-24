from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _write(dataframe, filename, format, **kwargs):
    """Write a PhyloPandas DataFrame to a sequence file.
    """
    seq_records = []
    # Iterate through dataframe rows
    for i, row in dataframe.iterrows():
        
        # Create a SeqRecord and append to list.
        seq_records.append(SeqRecord(
            Seq(row['sequence']),
            id=row['id'],
            name=row['name'],
            description=row['description']))
            
    # Write out sequences        
    SeqIO.write(seq_records, filename, format=format, **kwargs)
