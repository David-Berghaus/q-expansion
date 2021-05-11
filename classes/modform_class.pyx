"""
This class is just a placeholder and should be replaced by psage's "AutomorphicFormSpace" later
"""

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from classes.gamma_2_subgroup import Gamma_2_Subgroup

class ModForm(AutomorphicFormSpace):
    def __init__(self, G, weight=0):
        if isinstance(G, Gamma_2_Subgroup):
            self._group = G
            self._weight = weight
        else:
            super().__init__(G, weight)
    
    def weight(self):
        return self._weight

    def group(self):
        return self._group
