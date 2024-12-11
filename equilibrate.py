import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

# platform = Platform.getPlatformByName('CUDA')

dt = 0.002 * picoseconds   # Time step
tau_t = 1.0 * picoseconds  # Thermostat time constant (tautp)
cutoff = 10.0 * angstroms
report_interval = 5000      # Steps between output reports

# length of each stage
heat_time      = 3  * nanoseconds
nvt_time       = 2  * nanoseconds
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


def apply_posres_bb(weight, minim):
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

    if minim:
        positions = inpcrd.positions
    else:
        positions = simulation.state.getPositions()

    for atom, position in zip(prmtop.topology.atoms(), positions):
        if atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and atom.name in ['C','CA','N']:
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def apply_posres_sc(weight, minim):
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

    if minim:
        positions = inpcrd.positions
    else:
        positions = simulation.state.getPositions()

    for atom, position in zip(prmtop.topology.atoms(), positions):
        if atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and atom.element.symbol != "H" and atom.name not in ['C','CA','N']:
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def apply_posres_lig(weight, minim, lig_center_atoms):
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

    if minim:
        positions = inpcrd.positions
    else:
        positions = simulation.state.getPositions()

    for atom, position in zip(prmtop.topology.atoms(), positions):
        if atom.residue.name == 'LIG' and atom.element.symbol != "H" and atom.name not in lig_center_atoms:
            x0, y0, z0 = position
            restraint.addParticle(atom.index, [x0, y0, z0])

    # print(f"Number of restrained particles: {restraint.getNumParticles()}")
    force_index = system.addForce(restraint)

    # Reinitialize the simulation context to update forces
    simulation.context.reinitialize(preserveState=True)

    return(force_index, posres_name, weight)


def apply_posres_lig_center(weight, minim, lig_center_atoms):
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

    if minim:
        positions = inpcrd.positions
    else:
        positions = simulation.state.getPositions()

    for atom, position in zip(prmtop.topology.atoms(), positions):
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

    # adapted from https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html

    potential_energies = []
    constraint_energies = []
    constraint_strengths = []

    # you must override the report method and it must have this signature.
    def report(self, args):
        '''
        the report method is called every iteration of the minimization.
        '''

        potential_energy = args['system energy']
        constraint_energy = args['restraint energy']
        constraint_strength = args['restraint strength']

        self.potential_energies.append(potential_energy)
        self.constraint_strengths.append(constraint_strength)
        self.constraint_energies.append(constraint_energy)

        # The report method must return a bool specifying if minimization should be stopped.
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


def get_weight_list(posres_name, restraint_wt, start_time, stop_time, k0=5, decay_type='exponential'):
    """
    For a position restraint and the specified time points in the trajectory, make a list where len(list)==nsteps
    that determines the restraints weight at every step of the trajectory.
    Returns the weight list and also writes it to disk.
    """

    def exp_decay(a, j, k0, x):
        """
        This is the exp decay function, but normalized so that y reaches 0.
        the parameter k0 determines rate of decay, k is derived from it such that j does not influence it.
        i: the starting restraint weight
        j: x when y=0, corresponds to the timepoint when restraint should reach 0
        x: the xvalues, which corresponds to the simulation step numbers
        """
        import math

        k = k0/j
        y = a * ( math.exp(-k * x) - math.exp(-k * j) ) / ( 1 - math.exp(-k * j) )
        return(y)

    def linear_decay(a, j, x):
        """
        This will make just a line that starts at x=0, y=restr_weight
        and ends at x=end_time, y=0
        k0 is a placeholder
        """
        
        y = a * ( -x/j + 1)
        return(y)

    # stage 1
    stage1 = [ restraint_wt ] * int(np.round(start_time / dt))

    # stage 2
    stage2_numsteps = int(np.round( stop_time/dt - start_time/dt ))
    stage2_x = list(range(stage2_numsteps))

    line_params = [restraint_wt, stage2_numsteps]
    if decay_type == 'linear':
        decay_function = linear_decay
    elif decay_type == 'exponential':
        decay_function = exp_decay
        line_params.append(k0)

    stage2 = [decay_function(*line_params, x) for x in stage2_x]

    # stage 3
    stage3 = [0.0] * int( np.round(npt_restr_time/dt - stop_time/dt) )

    with open('restr_weights.csv', 'a+') as rstfile:
        rstfile.write('weight,name\n')
        rstfile.writelines(f'{str(i)},{posres_name}\n' for i in [*stage1, *stage2, *stage3])
    return([*stage1, *stage2, *stage3])


if __name__ == "__main__":

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
    posres_lig_center = apply_posres_lig_center(weight=10, lig_center_atoms=lig_center_atoms, minim=True)
    posres_bb = apply_posres_bb(weight=10, minim=True)
    posres_lig = apply_posres_lig(weight=10, lig_center_atoms=lig_center_atoms, minim=True)
    posres_sc = apply_posres_sc(weight=10, minim=True)

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
    posres_bb = apply_posres_bb(weight=10, minim=False)
    posres_sc = apply_posres_sc(weight=5, minim=False)
    posres_lig_center = apply_posres_lig_center(weight=10, lig_center_atoms=lig_center_atoms, minim=False)
    posres_lig = apply_posres_lig(weight=5, lig_center_atoms=lig_center_atoms, minim=False)

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

        wt_list = get_weight_list(posres_name, weight, start_time, stop_time, decay_type='exponential')
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