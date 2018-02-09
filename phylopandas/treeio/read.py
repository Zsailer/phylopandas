import pandas
import dendropy

def _doc_template(schema):
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


def _read(filename=None, data=None, schema=None, add_node_labels=True):
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
        tree = dendropy.Tree.get(path=filename, schema=schema)
    elif data is not None:
        tree = dendropy.Tree.get(data=data, schema=schema)
    else:
        raise Exception('No tree given?')

    # Maximum distance from root.
    tree.max_distance_from_root()

    # Initialize the data object.
    data = {'type':[],
            'id':[],
            'parent':[],
            'length':[],
            'label':[],
            'distance': []}

    # Add labels to internal nodes if set to true.
    if add_node_labels:
        for i, node in enumerate(tree.internal_nodes()):
            node.label = str(i)

    for node in tree.nodes():
        # Get node type
        if node.is_leaf():
            type_ = 'leaf'
            label = str(node.taxon.label)
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

    # Construct dataframe.
    df = pandas.DataFrame(data)
    return df


def read_newick(filename=None,
                data=None,
                add_node_labels=True):
    __doc__ = _doc_template('newick')
    return _read(filename=filename, data=data, schema='newick',
                 add_node_labels=True)


def read_nexus(filename=None,
               data=None,
               add_node_labels=True):
    __doc__ = _doc_template('nexus')
    return _read(filename=filename, data=data, schema='nexus',
                 add_node_labels=True)


def read_nexml(filename=None,
               data=None,
               add_node_labels=True):
    __doc__ = _doc_template('nexml')
    return _read(filename=filename, data=data, schema='nexml',
                 add_node_labels=True)
