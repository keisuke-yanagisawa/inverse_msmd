"""
Path utility functions.
"""
import os


def expandpath(path: str) -> str:
    """
    Expand ~ and $HOME and other environment variables.
    
    Parameters
    ----------
    path : str
        Path string that may contain ~ or environment variables
        
    Returns
    -------
    str
        Expanded path
    """
    path = os.path.expanduser(path)
    return os.path.expandvars(path)