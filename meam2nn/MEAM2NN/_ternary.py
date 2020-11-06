# coding: utf-8

import pandas as pd
import numpy as np

@property
def ternaries(self):
    """pandas.DataFrame: All MEAM ternary parameters."""
    return self.__ternaries

@property
def ternary_symbols(self):
    """list: All model symbols for which there are MEAM ternary paramters for."""
    symbollist = []
    for i in self.ternaries.index:
        ternary = self.ternaries.loc[i]
        symbollist.append([ternary.Sym1, ternary.Sym2, ternary.Sym3])
    return symbollist

def ternary(self, symbols):
    """
    The MEAM ternary parameters for a trio of symbol models.
    
    Parameters
    ----------
    symbols : str
        The symbols identifying the ternary model.  The symbols must be in the
        proper order.

    Returns
    -------
    pandas.Series
        The MEAM ternary parameters

    Raises
    ------
    ValueError
        If no parameters are found for the given symbols
    """
    try:
        return self.ternaries[((self.ternaries.Sym1==symbols[0]) & (self.ternaries.Sym2==symbols[1])  & (self.ternaries.Sym3==symbols[2]))].iloc[0]
    except:
        raise ValueError(f'MEAM parameters for ternary symbols {symbols} not found')

def load_ternaries(self, ternary_csv):
    """
    Loads MEAM ternary parameters from a csv and fills in fields if needed
    """
    # Load ternaries
    self.__ternaries = pd.read_csv(ternary_csv)

def lammps_parameter_ternary(self, ternary, i=1, j=2, k=3):
    """
    Generates the parameter file lines for a ternary interaction.

    Parameters
    ----------
    ternary : pd.Series
        The MEAM ternary parameters for a trio of symbol models.
    i : int, optional
        The index for the first symbol model: positive integer >= 1.  Default
        value is 1.
    j : int, optional
        The index for the second symbol model: positive integer > i.  Default
        value is 2.
    k : int, optional
        The index for the third symbol model: positive integer > j.  Default
        value is 3.

    Returns
    -------
    str
        The parameter lines
    """

    if not isinstance(i, int) or i <= 0:
        raise TypeError('i must be an integer > 0')
    if not isinstance(j, int) or j <= i:
        raise TypeError('j must be an integer > i')
    if not isinstance(k, int) or k <= j:
        raise TypeError('k must be an integer > j')

    contents =  f"Cmin({i},{j},{k}) = {ternary.Cmin_ikj:.15g}\n"
    #contents += f"Cmin({j},{i},{k}) = {ternary.Cmin_ikj}\n"
    contents += f"Cmin({i},{k},{j}) = {ternary.Cmin_ijk:.15g}\n"
    #contents += f"Cmin({k},{i},{j}) = {ternary.Cmin_ijk}\n"
    contents += f"Cmin({j},{k},{i}) = {ternary.Cmin_jik:.15g}\n"
    #contents += f"Cmin({k},{j},{i}) = {ternary.Cmin_jik}\n"

    contents += f"Cmax({i},{j},{k}) = {ternary.Cmax_ikj:.15g}\n"
    #contents += f"Cmax({j},{i},{k}) = {ternary.Cmax_ikj}\n"
    contents += f"Cmax({i},{k},{j}) = {ternary.Cmax_ijk:.15g}\n"
    #contents += f"Cmax({k},{i},{j}) = {ternary.Cmax_ijk}\n"
    contents += f"Cmax({j},{k},{i}) = {ternary.Cmax_jik:.15g}\n"
    #contents += f"Cmax({k},{j},{i}) = {ternary.Cmax_jik}\n"

    return contents