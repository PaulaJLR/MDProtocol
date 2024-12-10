from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds
from openmm import NoseHooverIntegrator
from config import SimulationConfig
from restraints import PositionRestraints


class Equilibration:

    def __init__(self, config:SimulationConfig):
        self.config = config
        
        self.prmtop = AmberPrmtopFile(config.top_name)
        self.inpcrd = AmberInpcrdFile(config.crd_name)
        self.system = self.prmtop.createSystem(
            nonbondedMethod=PME,
            nonbondedCutoff=config.cutoff,
            constraints=HBonds
        )
        self.integrator = NoseHooverIntegrator(
            config.start_temp,
            config.tau_t,
            config.dt
        )
        self.simulation = Simulation(
            self.prmtop.topology,
            self.system,
            self.integrator
        )
        self.simulation.context.setPositions(self.inpcrd.positions)

        self.position_restraints = []