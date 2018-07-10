import pandas
import dendropy


def _write(
    df,
    schema='newick',
    taxon_col=None,
    node_col=None,
    branch_lengths=True,
    ):
    """Write a phylopandas tree DataFrame to tree formats.
    """
    # Construct a list of nodes from dataframe.
    taxon_namespace = dendropy.TaxonNamespace()
    nodes = {}
    for idx in df.index:
        # Get node data.
        data = df.loc[idx]

        # Get taxon for node (if leaf node).
        taxon = None
        if data['type'] == 'leaf':
            if taxon_col is None:
                taxon = idx
            else:
                taxon = data[taxon_col]

            taxon = dendropy.Taxon(label=taxon)
            taxon_namespace.add_taxon(taxon)

        # Get label for node.
        if node_col is None:
            label = idx
        else:
            label = data[node_col]

        # Get edge length.
        edge_length = None
        if branch_lengths is True:
            edge_length = data['length']

        # Build a node
        n = dendropy.Node(
            taxon=taxon,
            label=label,
            edge_length=edge_length
        )

        nodes[idx] = n

    # Build branching pattern for nodes.
    root = None
    for idx, node in nodes.items():
        # Get node data.
        data = df.loc[idx]

        # Get children nodes
        children_idx = df[df['parent'] == data['id']].index
        children_nodes = [nodes[i] for i in children_idx]

        # Set child nodes
        nodes[idx].set_child_nodes(children_nodes)

        # Check if this is root.
        if data['parent'] is None:
            root = nodes[idx]

    # Build tree.
    tree = dendropy.Tree(
        seed_node=root,
        taxon_namespace=taxon_namespace
    )
    return tree
