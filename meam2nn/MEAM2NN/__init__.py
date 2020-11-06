# coding: utf-8

from ..tools import atomic_number, aslist

from importlib.resources import open_text
import pandas as pd

class MEAM2NN():
    """
    Class for handling and manipulating 2NN-MEAM parameters
    """
    from ._element import (elements, element_symbols, element, load_elements,
                           lammps_parameter_element)
    from ._alloy import (alloys, alloy_symbols, alloy, load_alloys,
                         lammps_parameter_alloy)
    from ._ternary import (ternaries, ternary_symbols, ternary, load_ternaries,
                           lammps_parameter_ternary)

    def __init__(self, element_csv=None, alloy_csv=None, ternary_csv=None):
        
        if element_csv is None:
            element_csv = open_text('meam2nn.parameters', 'element.csv')
        self.load_elements(element_csv)

        if alloy_csv is None:
            alloy_csv = open_text('meam2nn.parameters', 'alloy.csv')
        self.load_alloys(alloy_csv)

        if ternary_csv is None:
            ternary_csv = open_text('meam2nn.parameters', 'ternary.csv')
        self.load_ternaries(ternary_csv)

    def lammps_library(self, symbols):
        """
        Generate the MEAM library file contents

        Parameters
        ----------
        symbols : str or list, optional
            The element model symbols to include.  If not given, all
            element models will be included.
        """

        # Build the header
        contents =  '#elt\tlat\tz\tielement\tatwt\n'
        contents += '#alpha\t\tb0\tb1\tb2\tb3\talat\t\tesub\tasub\n'
        contents += '#t0\tt1\tt2\tt3\trozero\tibar\n'

        for symbol in aslist(symbols):
            element = self.element(symbol)
            contents += f"\n'{element.El}'\t'"
            contents += f"{element.lat}'\t"
            contents += f"{element.z}\t"
            contents += f"{element.ielement}\t"
            contents += f"{element.mass:.4f}\n"

            contents += f"{element.alpha:.10f}\t"
            contents += f"{element['beta(0)']:.3f}\t"
            contents += f"{element['beta(1)']:.3f}\t"
            contents += f"{element['beta(2)']:.3f}\t"
            contents += f"{element['beta(3)']:.3f}\t"
            contents += f"{element.alat:.10f}\t"
            contents += f"{element.Ec:.3f}\t"
            contents += f"{element.A:.3f}\n"

            contents += f"{element['t(0)']:.3f}\t"
            contents += f"{element['t(1)']:.3f}\t"
            contents += f"{element['t(2)']:.3f}\t"
            contents += f"{element['t(3)']:.3f}\t"
            contents += f"{element.Rho_zero:.3f}\t"
            contents += f"{element.ibar}\n"
        
        return contents
    
    def lammps_parameter(self, symbols, rcut=None, delr=0.1):
        """
        Generate the MEAM parameter file contents

        Parameters
        ----------
        symbols : str or list, optional
            The element model symbols to include.  If not given, all
            element models will be included.
        rcut : float, optional
            The cutoff distance to use for the model.  If not given, will use
            the largest set cutoff associated with the involved interactions.
        delr : float, optional
            The length of smoothing distance for cutoff function.  Default
            value is 0.1.
        """

        symbols = aslist(symbols)

        # Get involved interactions
        elements = []
        alloys = []
        alloys_ij = []
        ternaries = []
        ternaries_ijk = []

        for i in range(len(symbols)):
            isymbol = symbols[i]
            elements.append(self.element(isymbol))

            for j in range(i+1, len(symbols)):
                jsymbol = symbols[j]
                alloys.append(self.alloy([isymbol, jsymbol]))
                alloys_ij.append([i+1, j+1])

                for k in range(j+1, len(symbols)):
                    ksymbol = symbols[k]
                    ternaries.append(self.ternary([isymbol, jsymbol, ksymbol]))
                    ternaries_ijk.append([i+1, j+1, k+1])

        # Find default cutoff if needed
        if rcut is None:
            rcut = 0
            for element in elements:
                if element.Rcut > rcut:
                    rcut = element.Rcut
            for alloy in alloys:
                if alloy.Rcut > rcut:
                    rcut = alloy.Rcut
            for ternary in ternaries:
                if ternary.Rcut > rcut:
                    rcut = ternary.Rcut

        # Generate the global lines
        contents =  f"# Parameters for {' '.join(symbols)}\n"
        contents += f"rc = {rcut}\n"
        contents += f"delr = {delr}\n"
        contents += f"augt1 = 0\n"
        contents += f"erose_form = 2\n"
        contents += f"ialloy = 2\n"
        contents += f"mixture_ref_t = 0\n"
        contents += f"emb_lin_neg = 0\n"
        contents += f"bkgd_dyn = 0\n"

        # Generate the elemental lines
        for i, element in enumerate(elements):
            contents += '\n' + self.lammps_parameter_element(element, i+1)

        # Generate the alloy lines
        for alloy, ij in zip(alloys, alloys_ij):
            contents += '\n' + self.lammps_parameter_alloy(alloy, ij[0], ij[1])

        # Generate the ternary lines
        for ternary, ijk in zip(ternaries, ternaries_ijk):
            contents += '\n' + self.lammps_parameter_ternary(ternary, ijk[0], ijk[1], ijk[2])
        
        return contents