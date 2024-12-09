# config.py
from dataclasses import dataclass
from openmm.unit import *

@dataclass
class SimulationConfig:

    top_name:         str      = 'complex.prmtop'
    crd_name:         str      = 'complex.rst7'

    dt:               Quantity = 0.002 * picoseconds
    tau_t:            Quantity = 1.0 * picoseconds
    cutoff:           Quantity = 10.0 * angstroms
    report_interval:  int      = 5000
    
    heat_time:        Quantity = 3 * nanoseconds
    nvt_time:         Quantity = 2 * nanoseconds
    npt_restr_time:   Quantity = 15 * nanoseconds
    npt_time:         Quantity = 5 * nanoseconds
    
    start_temp:       Quantity = 100 * kelvin
    end_temp:         Quantity = 300 * kelvin
    
    pressure:         Quantity = 1 * bar

    lig_anchor_atoms:          = None

    def __post_init__(self):

        self.total_steps = (self.heat_time + self.nvt_time + self.npt_restr_time + self.npt_time) / self.dt


@dataclass
class RestraintConfig:

    minim_posres_bb_weight: Quantity = 10.0 * kilocalories_per_mole/angstroms**2
    minim_posres_sc_weight: Quantity = 10.0 * kilocalories_per_mole/angstroms**2
    minim_posres_lig_weight: Quantity = 10.0 * kilocalories_per_mole/angstroms**2
    minim_posres_lig_center_weight: Quantity = 10.0 * kilocalories_per_mole/angstroms**2

    posres_bb_weight: Quantity = 10.0 * kilocalories_per_mole/angstroms**2
    posres_sc_weight: Quantity = 5.0 * kilocalories_per_mole/angstroms**2
    posres_lig_center_weight: Quantity = 10.0 * kilocalories_per_mole/angstroms**2
    posres_lig_weight: Quantity = 5.0 * kilocalories_per_mole/angstroms**2

    lig_anchor_atoms = ['C1', 'C2', 'C3', 'C4', 'C5', 'O3', 'C9', 'C12', 'C13', 'C14', 'C15', 'O10', 'C20', 'C21', 'C22', 'C23', 'C24', 'O16', 'C28', 'C31', 'C32', 'C33', 'C34', 'O22']

    def posres_bb_mask(atom):
        return(
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and
            atom.name in ['C','CA','N']
        )

    def posres_sc_mask(atom):
        return(
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and
            atom.element.symbol != "H" and
            atom.name not in ['C','CA','N']
        )

    def posres_lig_anchor(atom, lig_anchor_atoms):
        return(
            atom.residue.name == 'LIG' and
            atom.name in lig_anchor_atoms
        )
    
    def posres_lig_ext(atom, lig_anchor_atoms):
        return(
            atom.residue.name == 'LIG' and
            atom.element.symbol != "H" and
            atom.name not in lig_anchor_atoms
        )