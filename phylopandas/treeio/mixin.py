

class DataFrameTreeMethods:
    """Tree operations.
    """

    def prune(self, **cell):
        """Prune a node or leaf from the current tree.
        
        Example call: 
            
            df.phylo.prune(id='1')
        
        Parameters
        ----------
        cell: 
            a column label and value pair. 
            
        Returns
        -------
        df : 
            DataFrame with the items in the subtree 
            below the given node pruned out.
        """
        if cell == {}:
            raise Exception("Only one option allowed.")
        if len(cell) > 1:
            raise Exception("Only one item allowed")
        
        df = self._data

        col = list(cell.keys())[0]
        val = list(cell.values())[0]
        
        # Convert to tree
        tree = df.phylo.to_dendropy(node_col=col)
        
        # Find node with the given label
        node = tree.find_node_with_label(val)
        
        # Get all descendant nodes
        descendants = list(node.postorder_iter())
        
        # Remove them from dataframe
        nodes_to_rm = [d.label for d in descendants]
        return df.loc[~df[col].isin(nodes_to_rm)]