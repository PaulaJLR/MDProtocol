from config import SimulationConfig, RestraintConfig
from simulation import Equilibration

import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

simconf = SimulationConfig(
    lig_center_atoms=['C1', 'C2', 'C3', 'C4', 'C5', 'O3', 'C9', 'C12', 'C13', 'C14', 'C15', 'O10', 'C20', 'C21', 'C22', 'C23', 'C24', 'O16', 'C28', 'C31', 'C32', 'C33', 'C34', 'O22']
)
restrconf = RestraintConfig()

simulation = Equilibration(simconf)

simulation.restraints.apply_restraints()

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


def apply_posres_lig(weight, lig_center_atoms):
    """
    Apply posres on the ligand, not center
    = :LIG and not H, and not lig_center_atoms
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
        if atom.residue.name == 'LIG' and atom.element.symbol != "H" and atom.name not in lig_center_atoms:
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def apply_posres_lig_center(weight, lig_center_atoms):
    """
    Apply posres on the ligand's center
    = :LIG and lig_center_atoms
    """

    posres_name = 'k_ligc'

    ## pos res heavy atoms solute
    force = weight * kilocalories_per_mole/angstroms**2
    restraint = CustomExternalForce(f'{posres_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint.addGlobalParameter(posres_name, force)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    for atom, position in zip(prmtop.topology.atoms(), inpcrd.positions):
        if atom.residue.name == 'LIG' and atom.name in lig_center_atoms:
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


# The class can have any name but it must subclass MinimizationReporter.
class MyMinimizationReporter(MinimizationReporter):

    # within the class you can declare variables that persist throughout the
    # minimization

    potential_energies = [] # array to record progress
    constraint_energies = []
    constraint_strengths = []

    # you must override the report method and it must have this signature.
    def report(self, iteration, x, grad, args):
        '''
        the report method is called every iteration of the minimization.

        Args:
            iteration (int): The index of the current iteration. This refers
                             to the current call to the L-BFGS optimizer.
                             Each time the minimizer increases the restraint strength,
                             the iteration index is reset to 0.

            x (array-like): The current particle positions in flattened order:
                            the three coordinates of the first particle,
                            then the three coordinates of the second particle, etc.

            grad (array-like): The current gradient of the objective function
                               (potential energy plus restraint energy) with
                               respect to the particle coordinates, in flattened order.

            args (dict): Additional statistics described above about the current state of minimization.
                         In particular:
                         “system energy”: the current potential energy of the system
                         “restraint energy”: the energy of the harmonic restraints
                         “restraint strength”: the force constant of the restraints (in kJ/mol/nm^2)
                         “max constraint error”: the maximum relative error in the length of any constraint

        Returns:
            bool : Specify if minimization should be stopped.
        '''

        # Within the report method you write the code you want to be executed at
        # each iteration of the minimization.
        # In this example we get the current energy, print it to the screen, and save it to an array.

        potential_energy = args['system energy']
        constraint_energy = args['restraint energy']
        constraint_strength = args['restraint strength']

        self.potential_energies.append(potential_energy)
        self.constraint_strengths.append(constraint_strength)
        self.constraint_energies.append(constraint_energy)

        # The report method must return a bool specifying if minimization should be stopped.
        # You can use this functionality for early termination.
        return False


def minimize(minim_name):

    reporter = MyMinimizationReporter()
    simulation.minimizeEnergy(reporter=reporter)
    save_rst7(top_name, f'{minim_name}.rst7')

    with open(f'{minim_name}.dat', 'a+') as minimf:
        minimf.write('potential_energy\tconstraint_energy\tconstraint_strength\n')
        for i in range(len(reporter.potential_energies)):
            potential_energy = float(reporter.potential_energies[i])
            constraint_energy = float(reporter.constraint_energies[i])
            constraint_strength = float(reporter.constraint_strengths[i])
            
            minimf.write(f'{potential_energy}\t{constraint_energy}\t{constraint_strength}\n')


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