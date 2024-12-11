from config import SimulationConfig, RestraintConfig
from simulation import Equilibration
from restraints import PositionRestraints

import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

simconf = SimulationConfig()
restrconf = RestraintConfig()

"""
prep system
"""
equilibration = Equilibration(simconf)

# add position restraints for minimization
## protein backbone
posres_bb = PositionRestraints(equilibration, 'k_bb', restrconf.minim_posres_bb_weight, restrconf.posres_bb_mask)
## ligand anchor
posres_anc = PositionRestraints(equilibration, 'k_la', restrconf.minim_posres_lig_anc_weight, restrconf.posres_lig_anc_mask)
## protein sidechain
posres_sc = PositionRestraints(equilibration, 'k_sc', restrconf.minim_posres_sc_weight, restrconf.posres_sc_mask)
## ligand extension
posres_ext = PositionRestraints(equilibration, 'k_le', restrconf.minim_posres_lig_ext_weight, restrconf.posres_lig_ext_mask)
# attention to order of adding restraints. each will have a force index,
# and then be removed by this index. If force a index is 0, and force b is 1,
# when i remove force a, force be will now be 0. This means that to make my life
# easier, the last one to be added has to be the first to go. if i am going to be
# removing them in the order a,b,c, then i need to add them in the order c,b,a.


"""
minimize
"""

# minimize 1st round
equilibration.minimize('minim1')

# remove posres sc:
posres_ext.remove()
posres_sc.remove()
# minimize 2nd round
equilibration.minimize('minim2')

# remove posres bb:
posres_anc.remove()
posres_bb.remove()
equilibration.minimize('minim3')



def add_reporters(out_name):

    simulation.reporters.clear()

    simulation.reporters.append(StateDataReporter(
        f"{out_name}.log", report_interval, step=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=total_steps, separator='\t'
    ))
    simulation.reporters.append(DCDReporter(f'{out_name}.dcd', report_interval))  # Trajectory output
    simulation.reporters.append(CheckpointReporter(f'{out_name}.chk', report_interval))  # Checkpointing


def heat():

    temps = np.linspace(start_temp._value, end_temp._value, nsteps)
    for i in range(len(temps)):
        temp = temps[i]
        integrator.setTemperature(temp * kelvin)
        simulation.step(1)


def add_barostat():

    pressure = 1 * bar

    # Add Monte Carlo barostat for pressure coupling
    barostat = MonteCarloBarostat(pressure, end_temp)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)


def get_weight_list(restraint_wt, start_time, stop_time):

    stage1 = [ restraint_wt ] * int(np.round(start_time / dt))
    stage2 = np.linspace( restraint_wt, 0.0, int(np.round( stop_time/dt - start_time/dt )) )
    stage3 = [0.0] * int( np.round(npt_restr_time/dt - stop_time/dt) )

    with open('restr_weights.txt', 'a+') as rstfile:
        rstfile.writelines(str(i)+'\n' for i in [*stage1, *stage2, *stage3])
        rstfile.write('===========================')
    return([*stage1, *stage2, *stage3])


"""
prep
"""

prep_system(top_name, crd_name)
lig_center_atoms = ['C1', 'C2', 'C3', 'C4', 'C5', 'O3', 'C9', 'C12', 'C13', 'C14', 'C15', 'O10', 'C20', 'C21', 'C22', 'C23', 'C24', 'O16', 'C28', 'C31', 'C32', 'C33', 'C34', 'O22']

"""
minimize
"""

# attention to order of adding restraints. each will have a force index,
# and then be removed by this index. If force a index is 0, and force b is 1,
# when i remove force a, force be will now be 0. This means that to make my life
# easier, the last one to be added has to be the first to go. if i am going to be
# removing them in the order a,b,c, then i need to add them in the order c,b,a.
posres_lig_center = apply_posres_lig_center(10, lig_center_atoms)
posres_bb = apply_posres_bb(10)
posres_lig = apply_posres_lig(10, lig_center_atoms)
posres_sc = apply_posres_sc(10)

# minimize
minimize('minim1')

# remove posres sc:
remove_posres(posres_sc[0])
remove_posres(posres_lig[0])
# minimize
minimize('minim2')

# remove posres bb:
remove_posres(posres_bb[0])
remove_posres(posres_lig_center[0])
# minimize
minimize('minim3')


"""
nvt
"""

### heating nvt
# Add reporters to output data
add_reporters('nvt_heat')

# posres protein bb at 10 and sc at 5
posres_bb = apply_posres_bb(10)
posres_sc = apply_posres_sc(5)
posres_lig_center = apply_posres_lig_center(10, lig_center_atoms)
posres_lig = apply_posres_lig(5, lig_center_atoms)

nsteps = int(heat_time / dt)
heat()

save_rst7(top_name, 'nvt_heat.rst7')


### nvt at temp
add_reporters('nvt')

nsteps = int(nvt_time / dt)
simulation.step(nsteps)

save_rst7(top_name, 'nvt.rst7')


"""
npt - reduce posres
"""

add_barostat()
add_reporters('npt_posres')

nsteps = int(np.round(npt_restr_time/dt))

#----- edit here -----#
# this uses the same weights from nvt
restraints = [
    
    #posres name, start time,  end time
    (posres_sc, 0*nanoseconds, npt_restr_time/3),
    (posres_bb, npt_restr_time/3, npt_restr_time),

    (posres_lig, 0*nanoseconds, npt_restr_time/3),
    (posres_lig_center, npt_restr_time/3, npt_restr_time),

]
#---------------------#

wt_lists = []

for restr in restraints:

    restraint, start_time, stop_time = restr
    force_index, posres_name, weight = restraint

    wt_list = get_weight_list(weight, start_time, stop_time)
    wt_lists.append(wt_list)


for i in range(nsteps):

    for j in range(len(restraints)):

        restraint, start_time, stop_time = restraints[j]
        force_index, posres_name, weight_0 = restraint
        weight = wt_lists[j][i]

        simulation.context.setParameter(posres_name, weight * kilocalories_per_mole/angstroms**2)

    simulation.step(1)

save_rst7(top_name, 'npt_posres.rst7')


"""
npt - no posres
"""

add_reporters('npt')

# no need to remove the restraints since they were brought to 0 in the prev step
nsteps = int(npt_time / dt)
simulation.step(nsteps)

save_rst7(top_name, 'npt1.rst7')