import pytest

from ..read import read_newick

# Newick strings pulled from https://en.wikipedia.org/wiki/Newick_format
NEWICKS = [
    "(,,(,));",                               # no nodes are named
    "(A,B,(C,D));",                           # leaf nodes are named
    "(A,B,(C,D)E)F;",                         # all nodes are named
    "(:0.1,:0.2,(:0.3,:0.4):0.5);",           # all but root node have a distance to parent
    "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",       # all have a distance to parent
    "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",       # distances and leaf names (popular)
    "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",     # distances and all names
    "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;"     # a tree rooted on a leaf node (rare)
]


def test_read_newick():
    for n in NEWICKS:
        print("Testing: {}".format(n))
        n = read_newick(data=n)

    assert True