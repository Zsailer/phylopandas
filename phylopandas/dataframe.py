
import pandas as pd
from .write import _write

class DataFrame(pd.DataFrame):
    
    def to_fasta(self, filename):
        """"""
        _write(self, filename, format="fasta")
        
    def to_phylip(self, filename):
        """"""
        _write(self, filename, format="phylip")
    
    
