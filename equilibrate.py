import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

platform = Platform.getPlatformByName('CUDA')

dt = 0.002 * picoseconds   # Time step
tau_t = 1.0 * picoseconds  # Thermostat time constant (tautp)
cutoff = 10.0 * angstroms
report_interval = 5000      # Steps between output reports

# length of each stage
heat_time      = 7  * nanoseconds
nvt_time       = 3  * nanoseconds
npt_restr_time = 15 * nanoseconds
npt_time       = 5  * nanoseconds

# temperature settings for heating
start_temp = 100 * kelvin
end_temp = 300 * kelvin

top_name = 'complex.prmtop'
crd_name = 'complex.rst7'

####
total_steps = (heat_time + nvt_time + npt_restr_time + npt_time) / dt
####


def prep_system(top_name, crd_name):

    global system
    global integrator
    global simulation
    global prmtop
    global inpcrd

    prmtop = AmberPrmtopFile(top_name)
    inpcrd = AmberInpcrdFile(crd_name)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff,
        constraints=HBonds  # SHAKE constraints on bonds involving hydrogen
    )

    integrator = NoseHooverIntegrator(start_temp, tau_t, dt)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)


def apply_posres_bb(weight):
    """
    Apply posres on protein backbone
    = (not water and not ions and not LIG) and (C,CA,N)
    """

    posres_name = 'k_bb'

    ## pos res heavy atoms solute
    force = weight * kilocalories_per_mole/angstroms**2
    restraint = CustomExternalForce(f'{posres_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint.addGlobalParameter(posres_name, force)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    for atom, position in zip(prmtop.topology.atoms(), inpcrd.positions):
        if atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and atom.name in ['C','CA','N']:
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def apply_posres_sc(weight):
    """
    Apply posres on protein side chain
    = (not water and not ions and not LIG) and not ('C','CA','N') and not H
    """

    posres_name = 'k_sc'

    force = weight * kilocalories_per_mole/angstroms**2
    restraint = CustomExternalForce(f'{posres_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint.addGlobalParameter(posres_name, force)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    for atom, position in zip(prmtop.topology.atoms(), inpcrd.positions):
        if atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and atom.element.symbol != "H" and atom.name not in ['C','CA','N']:
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def apply_posres_lig(weight):
    """
    Apply posres on the ligand
    = :LIG
    """

    posres_name = 'k_lig'

    ## pos res heavy atoms solute
    force = weight * kilocalories_per_mole/angstroms**2
    restraint = CustomExternalForce(f'{posres_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint.addGlobalParameter(posres_name, force)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    for atom, position in zip(prmtop.topology.atoms(), inpcrd.positions):
        if atom.residue.name == 'LIG':
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def remove_posres(force_index):

    # Remove the restraint force before continuing dynamics
    system.removeForce(force_index)
    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)


def save_rst7(top_name, out_crd_name):
    """
    use parmed to get positions and velocities and save to rst7
    """

    state = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    positions = state.getPositions()
    velocities = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()

    # from simtk.unit import nanometers, picoseconds, angstroms
    positions_in_angstroms = positions.value_in_unit(angstroms)
    velocities_in_angstroms_per_ps = velocities.value_in_unit(angstroms / picoseconds)

    amber_topology = pmd.load_file(top_name)
    amber_topology.positions = positions
    amber_topology.velocities = velocities
    amber_topology.save(out_crd_name, overwrite=True)


def minimize(minim_name):

    simulation.minimizeEnergy()
    save_rst7(top_name, f'{minim_name}.rst7')


def add_reporters(out_name):

    simulation.reporters.clear()

    simulation.reporters.append(StateDataReporter(
        f"{out_name}.log", report_interval, step=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=total_steps, separator='\t'
    ))
    simulation.reporters.append(DCDReporter(f'{out_name}.dcd', report_interval))  # Trajectory output
    simulation.reporters.append(CheckpointReporter(f'{out_name}.chk', report_interval))  # Checkpointing


def heat():

    temps = np.linspace(start_temp._value, end_temp._value, int(nsteps/temp_update))
    for i in range(len(temps)):
        temp = temps[i]
        integrator.setTemperature(temp * kelvin)
        simulation.step(temp_update)


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

"""
minimize
"""

# posres protein bb and sc:
posres_bb = apply_posres_bb(10)
posres_sc = apply_posres_sc(10)

# minimize
minimize('minim1')
save_rst7(top_name, 'minim1.rst7')

# remove posres sc:
remove_posres(posres_sc[0])
# minimize
minimize('minim2')
save_rst7(top_name, 'minim2.rst7')

# remove posres bb:
remove_posres(posres_bb[0])
# minimize
minimize('minim3')
save_rst7(top_name, 'minim3.rst7')


"""
nvt
"""

### heating nvt
# Add reporters to output data
add_reporters('nvt_heat')

# posres protein bb at 10 and sc at 5
posres_bb = apply_posres_bb(10)
posres_sc = apply_posres_sc(5)

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
    (posres_bb, npt_restr_time/3, npt_restr_time)
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

save_rst7(top_name, 'npt.rst7')