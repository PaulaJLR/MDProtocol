from dataclasses import dataclass, field
from openmm.unit import *
from typing import Callable
import math

@dataclass
class SimulationConfig:

    top_name:         str   = 'system.prmtop'
    crd_name:         str   = 'system.rst7'

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
        self.start_temp = self.start_temp * kelvin # has to be kelvin
        self.end_temp = self.end_temp * kelvin     # has to be kelvin
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
    lig_resname:      list     = field(default_factory=list)
    lig_anchor_atoms: list     = None
    structural_waters:list     = None
    mask_func_name:   str      = None
    decay_func:       Callable = None
    decay_rate:       float    = 5.0

    def __post_init__(self):

        self.weight = self.weight * kilocalories/mole/angstroms**2
        self.minim_weight = self.minim_weight * kilocalories/mole/angstroms**2
        self.start_time = self.start_time * nanoseconds
        self.end_time = self.end_time * nanoseconds

        if self.mask_func_name is None:
            self.mask_func_name = 'posres_bb_mask'
        if self.decay_func is None:
            self.decay_func = self.exp_decay
        
        mask_funcs = {
            'posres_bb_mask':self.posres_bb_mask,
            'posres_sc_mask':self.posres_sc_mask,
            'posres_liganc_mask':self.posres_liganc_mask,
            'posres_ligext_mask':self.posres_ligext_mask,
            'posres_water_mask':self.posres_water_mask
        }
        self.mask_func = mask_funcs[self.mask_func_name]


    def posres_bb_mask(self, atom):
        return(
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', *self.lig_resname] and
            atom.name in ['C','CA','N','O']
        )

    def posres_sc_mask(self, atom):
        return(
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', *self.lig_resname] and
            atom.element.symbol != "H" and
            atom.name not in ['C','CA','N','O']
        )

    def posres_liganc_mask(self, atom):
        return(
            atom.residue.name in self.lig_resname and
            atom.name in self.lig_anchor_atoms
        )
    
    def posres_ligext_mask(self, atom):
        return(
            atom.residue.name in self.lig_resname and
            atom.element.symbol != "H" and
            atom.name not in self.lig_anchor_atoms
        )
    
    def posres_water_mask(self, atom):
        return(
            (atom.residue.name == 'WAT' or atom.residue.name == 'HOH') and
            atom.residue.id in self.structural_waters and
            atom.name == 'O'
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
