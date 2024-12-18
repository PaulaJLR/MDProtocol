from config import SimulationConfig, RestraintConfig
from simulation import Equilibration
from restraints import PositionRestraints

from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd
from tools import add_barostat


"""
configure system
"""

# configure simulation
simconf = SimulationConfig(
    # add here all options that should not be kept default
    top_name='complex.prmtop',
    crd_name='complex.rst7',
)

lig_resname = 'LIG'
lig_anchor_atoms = ['C1', 'C2', 'C3', 'C4', 'C5', 'O3', 'C9', 'C12', 'C13', 'C14', 'C15', 'O10', 'C20', 'C21', 'C22', 'C23', 'C24', 'O16', 'C28', 'C31', 'C32', 'C33', 'C34', 'O22']
structural_waters = ['163', '164', '165', '166'] # residue numbers

# configure position restraints
config_posres_bb = RestraintConfig( # protein backbone
    name      = 'posres_bb',
    mask_func_name = 'posres_bb_mask'
)
config_posres_sc = RestraintConfig( # protein sidechain
    name       = 'posres_sc',
    weight     = 5.0,
    mask_func_name  = 'posres_sc_mask',
    start_time = 0,
    end_time   = simconf.get_value('npt_restr_time') / 3
)
config_posres_liganch = RestraintConfig( # ligand anchor
    name        = 'posres_liganc',
    mask_func_name   = 'posres_liganc_mask',
    lig_resname=lig_resname,
    lig_anchor_atoms=lig_anchor_atoms
)
config_posres_ligext = RestraintConfig( # ligand anchor
    name = 'posres_ligext',
    weight = 5.0,
    mask_func_name = 'posres_ligext_mask',
    lig_resname=lig_resname,
    lig_anchor_atoms=lig_anchor_atoms,
    start_time = 0,
    end_time = simconf.get_value('npt_restr_time') / 3
)
config_posres_wat = RestraintConfig(
    name = 'posres_wat',
    weight = 1.0,
    mask_func_name = 'posres_water_mask',
    structural_waters = structural_waters,
    start_time = 0,
    end_time = simconf.get_value('npt_restr_time') / 3
)

"""
prep system
"""

# start system
equilibration = Equilibration(simconf)

# initialize position restraints (this does not apply them yet)
posres_bb = PositionRestraints(equilibration, config_posres_bb)
posres_sc = PositionRestraints(equilibration, config_posres_sc)
posres_liganc = PositionRestraints(equilibration, config_posres_liganch)
posres_ligext = PositionRestraints(equilibration, config_posres_ligext)
posres_waters = PositionRestraints(equilibration, config_posres_wat)

"""
minimize
"""

# apply all restraints except struct waters
posres_bb.apply(minim=True)
posres_sc.apply(minim=True)
posres_liganc.apply(minim=True)
posres_ligext.apply(minim=True)
# minimize 1st round
equilibration.minimize('minim1')

# remove posres sc and ligext:
posres_sc.remove()
posres_ligext.remove()
# minimize 2nd round
equilibration.minimize('minim2')

# remove posres bb and liganc:
posres_bb.remove()
posres_liganc.remove()
# minimize 3rd round:
equilibration.minimize('minim3')


"""
nvt
"""

# start SHAKE
equilibration.apply_hbond_constraints()

# apply all position restraints for simulation
for restr in equilibration.position_restraints:
    restr.apply()

# increase temp linearly
equilibration.heat(simconf, 'nvt_heat')
# keep nvt at target temperatue
equilibration.nvt(simconf, 'nvt')


"""
npt - reduce posres
"""

add_barostat(simconf, equilibration)
equilibration.npt_posres(simconf, 'npt_posres')


"""
npt - no posres
"""

# even though restraint weights got to 0, I remove them
# to avoid unnecessery calculations:
for restr in equilibration.position_restraints:
    restr.remove()

equilibration.npt(simconf, 'npt')
