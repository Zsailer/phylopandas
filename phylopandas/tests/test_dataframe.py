import os
import pytest
from .. import dataframe
from . import path_to_dat, clean_dat

def test_read_fasta(path_to_dat):
    # Get path
    path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
    df = dataframe.read_fasta(path)
    
    # Tests
    keys = df.keys()
    assert type(df) == dataframe.DataFrame
    assert 'id' in keys
    assert 'sequence' in keys
    assert 'description' in keys
    assert 'name' in keys
    
def test_read_clustal(path_to_dat):
    # Get path
    path = os.path.join(path_to_dat, 'PF08793_seed.clustal')
    df = dataframe.read_clustal(path)
    
    # Tests
    keys = df.keys()
    assert type(df) == dataframe.DataFrame
    assert 'id' in keys
    assert 'sequence' in keys
    assert 'description' in keys
    assert 'name' in keys

def test_read_phylip(path_to_dat):
    # Get path
    path = os.path.join(path_to_dat, 'PF08793_seed.phylip')
    df = dataframe.read_phylip(path)
    
    # Tests
    keys = df.keys()
    assert type(df) == dataframe.DataFrame
    assert 'id' in keys
    assert 'sequence' in keys
    assert 'description' in keys
    assert 'name' in keys

class TestDataFrame(object):
    
    @pytest.mark.usefixtures("clean_dat")
    def test_to_csv(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = dataframe.read_fasta(path)
        
        # Write to csv
        csv_path = os.path.join(path_to_dat, 'test.csv')
        df.to_csv(csv_path)      
        assert os.path.exists(csv_path)
        
    @pytest.mark.usefixtures("clean_dat")
    def test_to_json(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = dataframe.read_fasta(path)
        
        # Write to csv
        json_path = os.path.join(path_to_dat, 'test.json')
        df.to_json(json_path)      
        assert os.path.exists(json_path)
        
    @pytest.mark.usefixtures("clean_dat")
    def test_to_fasta(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = dataframe.read_fasta(path)
        
        # Write to csv
        fasta_path = os.path.join(path_to_dat, 'test.fasta')
        df.to_fasta(fasta_path)      
        assert os.path.exists(fasta_path)
        
    @pytest.mark.usefixtures("clean_dat")
    def test_to_phylip(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = dataframe.read_fasta(path)
        
        # Write to csv
        phylip_path = os.path.join(path_to_dat, 'test.phylip')
        df.to_phylip(phylip_path)      
        assert os.path.exists(phylip_path)
        
    @pytest.mark.usefixtures("clean_dat")
    def test_to_embl(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = dataframe.read_fasta(path)
        
        # Write to csv
        embl_path = os.path.join(path_to_dat, 'test.embl')
        df.to_embl(alphabet='protein', filename=embl_path)      
        assert os.path.exists(embl_path)
        
    @pytest.mark.usefixtures("clean_dat")
    def test_to_nexus(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = dataframe.read_fasta(path)
        
        # Write to csv
        nexus_path = os.path.join(path_to_dat, 'test.nexus')
        df.to_nexus(alphabet='protein', filename=nexus_path)    
        assert os.path.exists(nexus_path)
