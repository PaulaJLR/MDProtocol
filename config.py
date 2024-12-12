from dataclasses import dataclass
from openmm.unit import *
from typing import Callable
import math

@dataclass
class SimulationConfig:

    top_name:         str      = 'system.prmtop'
    crd_name:         str      = 'system.rst7'

    dt:               float = 0.002
    tau_t:            float = 1.0
    cutoff:           float = 10.0
    report_interval:  int   = 5000
    
    heat_time:        float = 3.0
    nvt_time:         float = 2.0
    npt_restr_time:   float = 15.0
    npt_time:         float = 5.0
    
    start_temp:       float = 100.0
    end_temp:         float = 300.0
    
    pressure:         float = 1.0

    def __post_init__(self):

        # need to do this on post_init because cannot have mutable classes as default vals
        self.dt = self.dt * picoseconds
        self.tau_t = self.tau_t * picoseconds
        self.cutoff = self.cutoff * angstroms
        self.heat_time = self.heat_time * nanoseconds
        self.nvt_time = self.nvt_time * nanoseconds
        self.npt_restr_time = self.npt_restr_time * nanoseconds
        self.npt_time = self.npt_time * nanoseconds
        self.start_temp = self.start_temp * kelvin
        self.end_temp = self.end_temp * kelvin
        self.pressure = self.pressure * bar

        # get total steps
        self.total_steps = (self.heat_time + self.nvt_time + self.npt_restr_time + self.npt_time) / self.dt

    def get_value(self, config_item):
        """
        Returns the value of a config item which is an openmm Quantity in its appropriate unit
        """
        
        attribute = getattr(self, config_item)
        if type(attribute) == Quantity:
            return(attribute.value_in_unit(attribute.unit))
        else:
            return(attribute)


@dataclass
class RestraintConfig:

    weight:           float    = 10.0
    minim_weight:     float    = 10.0
    name:             str      = 'posres_bb'
    start_time:       float    = 5.0
    end_time:         float    = 15.0
    lig_anchor_atoms: list     = None
    mask_func:        Callable = None
    decay_func:       Callable = None
    decay_rate:       float    = 5.0

    def __post_init__(self):

        self.weight = self.weight * kilocalories/mole/angstroms**2
        self.minim_weight = self.minim_weight * kilocalories/mole/angstroms**2
        self.start_time = self.start_time * nanoseconds
        self.end_time = self.end_time * nanoseconds

        if self.mask_func is None:
            self.mask_func = self.posres_bb_mask
        if self.decay_func is None:
            self.decay_func = self.exp_decay

    @staticmethod
    def posres_bb_mask(atom):
        return(
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG', 'UNK'] and
            atom.name in ['C','CA','N']
        )

    @staticmethod
    def posres_sc_mask(atom):
        return(
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG', 'UNK'] and
            atom.element.symbol != "H" and
            atom.name not in ['C','CA','N']
        )

    @staticmethod
    def posres_liganc_mask(atom, lig_anchor_atoms):
        return(
            (atom.residue.name == 'LIG' or atom.residue.name == 'UNK') and
            atom.name in lig_anchor_atoms
        )
    
    @staticmethod
    def posres_ligext_mask(atom, lig_anchor_atoms):
        return(
            (atom.residue.name == 'LIG' or atom.residue.name == 'UNK') and
            atom.element.symbol != "H" and
            atom.name not in lig_anchor_atoms
        )
    

    @staticmethod
    def exp_decay(a, j, k0, x):
        """
        This is the exp decay function, but normalized so that y reaches 0.
        k0: determines rate of decay, k is derived from it such that j does not influence it.
        a: the starting restraint weight
        j: x when y=0, corresponds to the timepoint (num of steps) when restraint should reach 0
        x: the xvalues, which corresponds to the simulation step numbers
        """
        import math

        k = k0/j
        y = a * ( math.exp(-k * x) - math.exp(-k * j) ) / ( 1 - math.exp(-k * j) )
        return(y)
    

    @staticmethod
    def linear_decay(a, j, x):
        """
        This will make just a line that starts at x=0, y=restr_weight
        and ends at x=end_time, y=0
        """
        
        y = a * ( -x/j + 1)
        return(y)


    def get_value(self, config_item):
        """
        Returns the value of a config item which is an openmm Quantity in its appropriate unit
        """
        
        attribute = getattr(self, config_item)
        if type(attribute) == Quantity:
            return(attribute.value_in_unit(attribute.unit))
        else:
            return(attribute)