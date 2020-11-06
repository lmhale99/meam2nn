# coding: utf-8

import pandas as pd
import numpy as np

from ..tools import atomic_number

# Conversion constant between dyn/cm^2 to eV/A^3
b_scale = 0.624150907446076

@property
def elements(self):
    """pandas.DataFrame: All MEAM elemental parameters."""
    return self.__elements

@property
def element_symbols(self):
    """list: All model symbols for which there are MEAM elemental paramters for."""
    return self.elements.Sym.to_list()

def element(self, symbol):
    """
    The MEAM elemental parameters for a single symbol model.
    
    Parameters
    ----------
    symbol : str
        The symbol identifying the element model

    Returns
    -------
    pandas.Series
        The MEAM elemental parameters

    Raises
    ------
    ValueError
        If no parameters are found for the given symbol model
    """
    try:
        return self.elements[self.elements.Sym == symbol].iloc[0]
    except:
        raise ValueError(f'MEAM parameters for element symbol {symbol} not found')

def load_elements(self, element_csv):
    """
    Loads MEAM elemental parameters from a csv and fills in fields if needed
    """
    # Load elements
    self.__elements = pd.read_csv(element_csv)

    # Specify dicts of 1:1 values to assign
    ref_lat = {'FCC_A1':'fcc', 'BCC_A2':'bcc', 'DIA_A4':'dia', 'HCP_A3':'hcp', 'DIMER':'dim'}
    ref_z = {'FCC_A1':12, 'BCC_A2':8, 'DIA_A4':4, 'HCP_A3':12, 'DIMER':1}

    # initialize empty lists for values
    lats = []
    zs = []
    ielements = []
    t0s = []
    ibars = []
    alats = []
    bs = [] # Bulk mod in eV/A^3
    bstars = [] # Bulk mod in dyn/cm^2
    alphas = []
    omegas = []

    for i in self.elements.index:
        element = self.elements.loc[i]

        # Set 1:1 values
        lats.append(ref_lat[element['Ref.St']])
        zs.append(ref_z[element['Ref.St']])
        ielements.append(atomic_number(element.El))
        t0s.append(1.0)
        ibars.append(3)

        # Set dimer parameters (no real ref crystal)
        if element['Ref.St'] == 'DIMER':
            alat = element.Re
            if 'B' in element and pd.notna(element.B):
                b = element.B
                bstar = b / b_scale
                alpha = b
            elif 'B*' in element and pd.notna(element['B*']):
                bstar = element['B*']
                b = bstar * b_scale
                alpha = b
            elif 'alpha' in element and pd.notna(element.alpha):
                alpha = element.alpha
                b = alpha
                bstar = b / b_scale
            else:
                raise ValueError('B, B* or alpha must be specified')
            omega = element.Ec * b / 9

        else:
            # Set crystal structure parameters
            if element['Ref.St'] == 'FCC_A1':
                alat = element.Re * 2**0.5
                omega = alat**3 / 4
            elif element['Ref.St'] == 'BCC_A2':
                alat = element.Re * 2 * 3**0.5 / 3
                omega = alat**3 / 2 
            elif element['Ref.St'] == 'DIA_A4':
                alat = element.Re / (3**0.5 / 4)
                omega = alat**3 / 8
            elif element['Ref.St'] == 'HCP_A3':
                alat = element.Re
                omega = alat**3 / 2**0.5

            # Set bulk modulus/alpha parameters
            if 'B' in element and pd.notna(element.B):
                b = element.B
                bstar = b / b_scale
                alpha = (9 * b * omega / element.Ec)**0.5
            elif 'B*' in element and pd.notna(element['B*']):
                bstar = element['B*']
                b = bstar * b_scale
                alpha = (9 * b * omega / element.Ec)**0.5
            elif 'alpha' in element and pd.notna(element.alpha):
                alpha = element.alpha
                b = alpha**2 * element.Ec / (9 * omega)
                bstar = b / b_scale
            else:
                raise ValueError('B, B* or alpha must be specified')

        # Append values
        alats.append(alat)
        bs.append(b)
        bstars.append(bstar)
        alphas.append(alpha)
        omegas.append(omega)

    # Save values to DataFrame
    self.elements['lat'] = lats
    self.elements['z'] = zs
    self.elements['ielement'] = ielements
    self.elements['t(0)'] = t0s
    self.elements['ibar'] = ibars
    self.elements['alat'] = alats
    self.elements['B'] = bs
    self.elements['B*'] = bstars
    self.elements['alpha'] = alphas
    self.elements['omega'] = omegas

def lammps_parameter_element(self, element, i=1):
    """
    Generates the parameter file lines for an element.

    Parameters
    ----------
    element : pd.Series
        The MEAM elemental parameters for a single symbol model.
    i : int, optional
        The index for the symbol model: positive integer >= 1.  Default value
        is 1.

    Returns
    -------
    str
        The parameter lines
    """
    if not isinstance(i, int) or i <=0:
        raise TypeError('i must be an integer > 0')

    contents =  f"zbl({i},{i}) = 0\n"
    contents += f"nn2({i},{i}) = {element.nn2}\n"
    contents += f"rho0({i}) = {element.Rho_zero}\n"
    contents += f"Ec({i},{i}) = {element.Ec:.8g}\n"
    contents += f"re({i},{i}) = {element.Re:.8g}\n"
    contents += f"repuls({i},{i}) = {element['d-']:.6g}\n"
    contents += f"attrac({i},{i}) = {element['d+']:.6g}\n"
    contents += f"Cmin({i},{i},{i}) = {element.Cmin:.15g}\n"
    contents += f"Cmax({i},{i},{i}) = {element.Cmax:.15g}\n"

    return contents