__doc__ = """\
This module coverts the DataFrames returned by all of Pandas `read_` methods 
to PhyloPandas DataFrames. 
"""
import pandas
from functools import wraps
from .dataframe import DataFrame

def _pandas_to_phylopandas_converter(func):
    """Decorate to covert all Pandas `read_` functions to output phylopandas.DataFrame"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        df = func(*args, **kwargs)
        return DataFrame(df)
    return wrapper

def _convert_pandas_functions_to_phylopandas(pandas_module):
    """Convert all pandas functions to phylopandas functions"""
    for name in dir(pandas_module):
        # Iterate through pandas functions   
        if name[:5] == 'read_':
            func = getattr(pandas_module, name)
            setattr(pandas_module, name, _pandas_to_phylopandas_converter(func))

# Do the conversion
_convert_pandas_functions_to_phylopandas(pandas)
