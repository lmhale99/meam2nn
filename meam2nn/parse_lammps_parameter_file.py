from pathlib import Path

def parse_lammps_parameter_file(paramfile):
    """
    Parses all values from a LAMMPS MEAM parameter file.

    Parameters
    ----------
    paramfile : str, Path, or file-like object
        The parameter file name, contents or open file to read and parse/

    Returns
    -------
    dict
        The parsed MEAM parameters.
    """

    lines = None
    needtoclose = False

    # Check if str paramfile is a file name or file contents
    if isinstance(paramfile, str):
        try:
            assert Path(paramfile).is_file()
        except:
            lines = paramfile.split('\n')
        else:
            paramfile = Path(paramfile)
    
    # Open file if paramfile is a path
    if isinstance(paramfile, Path):
        paramfile = open(paramfile)
        needtoclose = True
    
    # Read file
    if hasattr(paramfile, 'readlines'):
        lines = paramfile.readlines()
        if needtoclose:
            paramfile.close()
    if lines is None:
        raise TypeError('paramfile not recognized as str contents, file path or file-like object')
    
    # Parse contents line by line
    params = {}
    for line in lines:
        line = line.strip()
        
        # Ignore empty and comment lines
        if len(line) == 0 or line[0] == '#':
            continue

        # Split key and value
        key, value = line.split('=')
        key = key.strip()
        value = value.strip()
        
        # Convert to int and float if possible
        try:
            value = int(value)
        except:
            try:
                value = float(value)
            except:
                pass
        params[key] = value
        
    return params