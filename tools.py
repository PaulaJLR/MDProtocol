import parmed as pmd

def save_rst7(equilibration, out_crd_name):
    """
    use parmed to get positions and velocities and save to rst7
    """

    state = equilibration.simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    positions = state.getPositions()
    velocities = state.getVelocities()

    amber_topology = pmd.openmm.load_topology(equilibration.prmtop.topology, equilibration.system)
    amber_topology.positions = positions
    amber_topology.velocities = velocities
    amber_topology.save(out_crd_name, overwrite=True)
