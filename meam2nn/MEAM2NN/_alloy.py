# coding: utf-8

import pandas as pd
import numpy as np

# Conversion constant between dyn/cm^2 to eV/A^3
b_scale = 0.624150907446076

@property
def alloys(self):
    """pandas.DataFrame: All MEAM alloy parameters."""
    return self.__alloys

@property
def alloy_symbols(self):
    """list: All model symbols for which there are MEAM alloy paramters for."""
    symbollist = []
    for i in self.alloys.index:
        alloy = self.alloys.loc[i]
        symbollist.append([alloy.Sym1, alloy.Sym2])
    return symbollist

def alloy(self, symbols):
    """
    The MEAM alloy parameters for a pair of symbol models.
    
    Parameters
    ----------
    symbols : str
        The symbols identifying the alloy model.  The symbols must be in the
        proper order.

    Returns
    -------
    pandas.Series
        The MEAM alloy parameters

    Raises
    ------
    ValueError
        If no parameters are found for the given symbols
    """
    try:
        return self.alloys[((self.alloys.Sym1==symbols[0]) & (self.alloys.Sym2==symbols[1]))].iloc[0]
    except:
        raise ValueError(f'MEAM parameters for alloy symbols {symbols} not found')

def load_alloys(self, alloy_csv):
    """
    Loads MEAM alloy parameters from a csv and fills in fields if needed
    """
    # Load alloys
    self.__alloys = pd.read_csv(alloy_csv)

    # Specify dicts of 1:1 values to assign
    ref_lat = {'FCC_B1':'b1', 'BCC_B2':'b2', 'ZnS_B3':'dia', 'L12A3B':'l12'}

    # initialize empty lists for values
    lats = []
    energies = []
    alats = []
    omegas = []
    bs = [] # Bulk mod in eV/A^3
    bstars = [] # Bulk mod in dyn/cm^2
    alphas = []

    for a in self.alloys.index:
        
        # Fetch alloy and associated elements
        alloy = self.alloys.loc[a]
        ielement = self.element(alloy.Sym1)
        jelement = self.element(alloy.Sym2)

        # Set 1:1 values
        lats.append(ref_lat[alloy['Ref.St']])

        # Set lattice-specific averaging fractions
        if alloy['Ref.St'] == 'L12A3B':
            iscale = 0.75
            jscale = 0.25
        else:
            iscale = jscale = 0.5

        # Compute Ec based on elements and delta_Ec
        Ec = iscale * ielement.Ec + jscale * jelement.Ec - alloy.delta_Ec
        energies.append(Ec)

        # Set d values using weighted average of elements if not given
        if pd.isna(alloy['d+']):
            self.alloys.loc[a, 'd+'] = iscale * ielement['d+'] + jscale * jelement['d+']
        if pd.isna(alloy['d-']):
            self.alloys.loc[a, 'd-'] = iscale * ielement['d-'] + jscale * jelement['d-']

        # Set crystal structure parameters
        re = alloy.Re

        # Use weighted average omega if re not given
        if pd.isna(re):
            omega = iscale * ielement.omega + jscale * jelement.omega
            if alloy['Ref.St'] == 'FCC_B1':
                alat = (8 * omega) ** (1/3)
                re = alat / 2
            elif alloy['Ref.St'] == 'BCC_B2':
                alat = (2 * omega)**(1/3)
                re = alat * 3**0.5 / 2
            elif alloy['Ref.St'] == 'ZnS_B3':
                alat = (8 * omega)**(1/3)
                re = alat * 3**0.5 / 4
            elif alloy['Ref.St'] == 'L12A3B':
                alat = (4 * omega)**(1/3)
                re = alat / 2**0.5
            self.alloys.loc[a, 'Re'] = re

        # Compute alat and omega using re and structure
        else:
            if alloy['Ref.St'] == 'FCC_B1':
                alat = 2 * re
                omega = alat ** 3 / 8
            elif alloy['Ref.St'] == 'BCC_B2':
                alat = re * 2 * 3**0.5 / 3
                omega = alat**3 / 2
            elif alloy['Ref.St'] == 'ZnS_B3':
                alat = 4 * re / 3**0.5
                omega = alat**3 / 8
            elif alloy['Ref.St'] == 'L12A3B':
                alat = 2**0.5 * re
                omega = alat**3 / 4
        alats.append(alat)
        omegas.append(omega)

        # Set bulk modulus/alpha parameters
        if 'B' in alloy and pd.notna(alloy.B):
            b = alloy.B
            bstar = b / b_scale
            alpha = (9 * b * omega / Ec)**0.5
        elif 'B*' in alloy and pd.notna(alloy['B*']):
            bstar = alloy['B*']
            b = bstar * b_scale
            alpha = (9 * b * omega / Ec)**0.5
        elif 'alpha' in alloy and pd.notna(alloy.alpha):
            alpha = alloy.alpha
            b = alpha**2 * Ec / (9 * omega)
            bstar = b / b_scale
        else:
            # Use weighted average if none given
            b = iscale * ielement.B + jscale * jelement.B
            bstar = b / b_scale
            alpha = (9 * b * omega / Ec)**0.5
        bs.append(b)
        bstars.append(bstar)
        alphas.append(alpha)

        # Set default Cmin values
        if pd.isna(alloy.Cmin_iji) or alloy.Cmin_iji == 0:
            self.alloys.loc[a, 'Cmin_iji'] = ielement.Cmin
        if pd.isna(alloy.Cmin_jij)or alloy.Cmin_jij == 0:
            self.alloys.loc[a, 'Cmin_jij'] = jelement.Cmin
        if pd.isna(alloy.Cmin_iij)or alloy.Cmin_iij == 0:
            self.alloys.loc[a, 'Cmin_iij'] = (0.5 * ielement.Cmin**0.5 + 0.5 * jelement.Cmin**0.5)**2
        if pd.isna(alloy.Cmin_ijj)or alloy.Cmin_ijj == 0:
            self.alloys.loc[a, 'Cmin_ijj'] = (0.5 * ielement.Cmin**0.5 + 0.5 * jelement.Cmin**0.5)**2

        # Set default Cmax values
        if alloy.Cmax_iji == 0:
            self.alloys.loc[a, 'Cmax_iji'] = 2.8
        elif pd.isna(alloy.Cmax_iji):
            self.alloys.loc[a, 'Cmax_iji'] = ielement.Cmax
        if alloy.Cmax_jij == 0:
            self.alloys.loc[a, 'Cmax_jij'] = 2.8
        elif pd.isna(alloy.Cmax_jij):
            self.alloys.loc[a, 'Cmax_jij'] = jelement.Cmax
        if alloy.Cmax_iij == 0:
            self.alloys.loc[a, 'Cmax_iij'] = 2.8
        elif pd.isna(alloy.Cmax_iij):
            self.alloys.loc[a, 'Cmax_iij'] = (0.5 * ielement.Cmax**0.5 + 0.5 * jelement.Cmax**0.5)**2
        if alloy.Cmax_ijj == 0:
            self.alloys.loc[a, 'Cmax_ijj'] = 2.8
        elif pd.isna(alloy.Cmax_ijj):
            self.alloys.loc[a, 'Cmax_ijj'] = (0.5 * ielement.Cmax**0.5 + 0.5 * jelement.Cmax**0.5)**2
    
    # Save values to DataFrame
    self.alloys['lat'] = lats
    self.alloys['Ec'] = energies
    self.alloys['alat'] = alats
    self.alloys['omega'] = omegas
    self.alloys['B'] = bs
    self.alloys['B*'] = bstars
    self.alloys['alpha'] = alphas

def lammps_parameter_alloy(self, alloy, i=1, j=2):
    """
    Generates the parameter file lines for an alloy interaction.

    Parameters
    ----------
    alloy : pd.Series
        The MEAM alloy parameters for a pair of symbol models.
    i : int, optional
        The index for the first symbol model: positive integer >= 1.  Default
        value is 1.
    j : int, optional
        The index for the second symbol model: positive integer > i.  Default
        value is 2.

    Returns
    -------
    str
        The parameter lines
    """

    if not isinstance(i, int) or i <= 0:
        raise TypeError('i must be an integer > 0')
    if not isinstance(j, int) or j <= i:
        raise TypeError('j must be an integer > i')

    contents =  f"zbl({i},{j}) = 0\n"
    contents += f"nn2({i},{j}) = 1\n"
    contents += f"Ec({i},{j}) = {alloy.Ec:.8g}\n"
    contents += f"re({i},{j}) = {alloy.Re:.8g}\n"
    contents += f"alpha({i},{j}) = {alloy.alpha:.15g}\n"
    contents += f"attrac({i},{j}) = {alloy['d+']:.6g}\n"
    contents += f"repuls({i},{j}) = {alloy['d-']:.6g}\n"

    contents += f"Cmin({i},{i},{j}) = {alloy.Cmin_iji:.15g}\n"
    contents += f"Cmin({j},{j},{i}) = {alloy.Cmin_jij:.15g}\n"
    contents += f"Cmin({i},{j},{i}) = {alloy.Cmin_iij:.15g}\n"
    contents += f"Cmin({i},{j},{j}) = {alloy.Cmin_ijj:.15g}\n"
    contents += f"Cmin({j},{i},{i}) = {alloy.Cmin_iij:.15g}\n"
    contents += f"Cmin({j},{i},{j}) = {alloy.Cmin_ijj:.15g}\n"

    contents += f"Cmax({i},{i},{j}) = {alloy.Cmax_iji:.15g}\n"
    contents += f"Cmax({j},{j},{i}) = {alloy.Cmax_jij:.15g}\n"
    contents += f"Cmax({i},{j},{i}) = {alloy.Cmax_iij:.15g}\n"
    contents += f"Cmax({i},{j},{j}) = {alloy.Cmax_ijj:.15g}\n"
    contents += f"Cmax({j},{i},{i}) = {alloy.Cmax_iij:.15g}\n"
    contents += f"Cmax({j},{i},{j}) = {alloy.Cmax_ijj:.15g}\n"

    contents += f"lattce({i},{j}) = '{alloy.lat}'\n"

    return contents