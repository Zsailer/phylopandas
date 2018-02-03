import os
import pytest
import phylopandas as ph
from . import path_to_dat, clean_dat


def test_read_fasta(path_to_dat):
    # Get path
    path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
    df = ph.read_fasta(path)

    # Tests
    keys = df.keys()
    assert type(df) == ph.DataFrame
    assert 'id' in keys
    assert 'sequence' in keys
    assert 'description' in keys
    assert 'name' in keys


def test_read_clustal(path_to_dat):
    # Get path
    path = os.path.join(path_to_dat, 'PF08793_seed.clustal')
    df = ph.read_clustal(path)

    # Tests
    keys = df.keys()
    assert type(df) == ph.DataFrame
    assert 'id' in keys
    assert 'sequence' in keys
    assert 'description' in keys
    assert 'name' in keys


def test_read_phylip(path_to_dat):
    # Get path
    path = os.path.join(path_to_dat, 'PF08793_seed.phylip')
    df = ph.read_phylip(path)

    # Tests
    keys = df.keys()
    assert type(df) == ph.DataFrame
    assert 'id' in keys
    assert 'sequence' in keys
    assert 'description' in keys
    assert 'name' in keys


class Testframe(object):

    @pytest.mark.usefixtures("clean_dat")
    def test_to_fasta(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = ph.read_fasta(path)
        print(df.phylo)
        # Write to fasta
        fasta_path = os.path.join(path_to_dat, 'test.fasta')
        df.phylo.to_fasta(fasta_path)
        assert os.path.exists(fasta_path)

    @pytest.mark.usefixtures("clean_dat")
    def test_to_phylip(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = ph.read_fasta(path)

        # Write to csv
        phylip_path = os.path.join(path_to_dat, 'test.phylip')
        df.phylo.to_phylip(phylip_path)
        assert os.path.exists(phylip_path)

    @pytest.mark.usefixtures("clean_dat")
    def test_to_embl(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = ph.read_fasta(path)

        # Write to csv
        embl_path = os.path.join(path_to_dat, 'test.embl')
        df.phylo.to_embl(alphabet='protein', filename=embl_path)
        assert os.path.exists(embl_path)

    @pytest.mark.usefixtures("clean_dat")
    def test_to_nexus(self, path_to_dat):
        path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
        df = ph.read_fasta(path)

        # Write to csv
        nexus_path = os.path.join(path_to_dat, 'test.nexus')
        df.phylo.to_nexus(alphabet='protein', filename=nexus_path)
        assert os.path.exists(nexus_path)
