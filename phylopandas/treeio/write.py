import pandas
import dendropy


def _write(df, schema):
    """Write a phylopandas tree DataFrame to tree format.s
    """
    # Get taxon namespace
    leaf_nodes = df[df["type"] == "leaf"]
    taxon_namespace = dendropy.TaxonNamespace(leaf_nodes["name"])

    # Initialize a tree
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)

    # Build tree
    node = {}
