import pandas
import dendropy
from ..utils import get_random_id

def _read_doc_template(schema):
    doc = """
    Read a {} tree into a phylopandas.DataFrame.

    The resulting DataFrame has the following columns:
        - name: label for each taxa or node.
        - id: unique id (created by phylopandas) given to each node.
        - type: type of node (leaf, internal, or root).
        - parent: parent id. necessary for constructing trees.
        - length: length of branch from parent to node.
        - distance: distance from root.

    Parameters
    ----------
    filename: str (default is None)
        {} file to read into DataFrame.

    data: str (default is None)
        {} string to parse and read into DataFrame.

    add_node_labels: bool
        If true, labels the internal nodes with numbers.

    Returns
    -------
    df: phylopandas.DataFrame
    """.format(schema, schema, schema)
    return doc


def _dendropy_to_dataframe(
    tree,
    add_node_labels=True,
    use_uids=True):
    """Convert Dendropy tree to Pandas dataframe."""
    # Maximum distance from root.
    tree.max_distance_from_root()

    # Initialize the data object.
    idx = []
    data = {
        'type': [],
        'id': [],
        'parent': [],
        'length': [],
        'label': [],
        'distance': []}

    if use_uids:
        data['uid'] = []

    # Add labels to internal nodes if set to true.
    if add_node_labels:
        for i, node in enumerate(tree.internal_nodes()):
            node.label = str(i)

    for node in tree.nodes():
        # Get node type
        if node.is_leaf():
            type_ = 'leaf'
            label = str(node.taxon.label).replace(' ', '_')
        elif node.is_internal():
            type_ = 'node'
            label = str(node.label)

        # Set node label and parent.
        id_ = label
        parent_node = node.parent_node
        length = node.edge_length
        distance = node.distance_from_root()

        # Is this node a root?
        if parent_node is None and length is None:
            parent_label = None
            parent_node = None
            length = 0
            distance = 0
            type_ = 'root'

        # Set parent node label
        elif parent_node.is_internal():
            parent_label = str(parent_node.label)

        else:
            raise Exception("Subtree is not attached to tree?")

        # Add this node to the data.
        data['type'].append(type_)
        data['id'].append(id_)
        data['parent'].append(parent_label)
        data['length'].append(length)
        data['label'].append(label)
        data['distance'].append(distance)

        if use_uids:
            data['uid'].append(get_random_id(10))

    # Construct dataframe.
    df = pandas.DataFrame(data)
    return df


def _read(
    filename=None,
    data=None,
    schema=None,
    add_node_labels=True,
    use_uids=True
    ):
    """Read a phylogenetic tree into a phylopandas.DataFrame.

    The resulting DataFrame has the following columns:
        - name: label for each taxa or node.
        - id: unique id (created by phylopandas) given to each node.
        - type: type of node (leaf, internal, or root).
        - parent: parent id. necessary for constructing trees.
        - length: length of branch from parent to node.
        - distance: distance from root.

    Parameters
    ----------
    filename: str (default is None)
        newick file to read into DataFrame.

    data: str (default is None)
        newick string to parse and read into DataFrame.

    add_node_labels: bool
        If true, labels the internal nodes with numbers.

    Returns
    -------
    df: phylopandas.DataFrame.
    """
    if filename is not None:
        # Use Dendropy to parse tree.
        tree = dendropy.Tree.get(
            path=filename,
            schema=schema,
            preserve_underscores=True)
    elif data is not None:
        tree = dendropy.Tree.get(
            data=data,
            schema=schema,
            preserve_underscores=True)
    else:
        raise Exception('No tree given?')

    df = _dendropy_to_dataframe(
        tree, 
        add_node_labels=add_node_labels,
        use_uids=use_uids
    )
    return df


def _read_method(schema):
    """Add a write method for named schema to a class.
    """
    def func(
        self,
        filename=None,
        data=None,
        add_node_labels=True,
        combine_on='index',
        use_uids=True,
        **kwargs):
        # Use generic write class to write data.
        df0 = self._data
        df1 = _read(
            filename=filename,
            data=data,
            schema=schema,
            add_node_labels=add_node_labels,
            use_uids=use_uids,
            **kwargs
        )
        return df0.phylo.combine(df1, on=combine_on)

    # Update docs
    func.__doc__ = _read_doc_template(schema)
    return func


def _read_function(schema):
    """Add a write method for named schema to a class.
    """
    def func(
        filename=None,
        data=None,
        add_node_labels=True,
        use_uids=True,
        **kwargs):
        # Use generic write class to write data.
        return _read(
            filename=filename,
            data=data,
            schema=schema,
            add_node_labels=add_node_labels,
            use_uids=use_uids,
            **kwargs
        )
    # Update docs
    func.__doc__ = _read_doc_template(schema)
    return func


def read_dendropy(
    df,         
    add_node_labels=True,
    use_uids=True):
    __doc__ = _read_doc_template('dendropy')

    df = _dendropy_to_dataframe(
        tree,
        add_node_labels=add_node_labels,
        use_uids=use_uids
    )
    return df

read_newick = _read_function('newick')
read_nexml = _read_function('nexml')
read_nexus_tree = _read_function('nexus')
