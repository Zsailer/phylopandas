import os
import pytest

import phylopandas as ph
from . import path_to_dat, clean_dat  # noqa


@pytest.mark.usefixtures("clean_dat")
def test_to_fasta(path_to_dat):
    path = os.path.join(path_to_dat, 'PF08793_seed.fasta')
    df = ph.read_fasta(path)

    # Extract a single row
    row = df.iloc[0]

    # Write row to fasta
    fasta_path = os.path.join(path_to_dat, 'test.fasta')
    row.phylo.to_fasta(fasta_path)
    assert os.path.exists(fasta_path)
