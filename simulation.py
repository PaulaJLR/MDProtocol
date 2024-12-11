from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds
from openmm import NoseHooverIntegrator
from config import SimulationConfig
import parmed as pmd
from tools import save_rst7
from reporters import MyMinimizationReporter

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

    def minimize(self, minim_name):

        reporter = MyMinimizationReporter()
        self.simulation.minimizeEnergy(reporter=reporter)
        save_rst7(self, f'{minim_name}.rst7')

        with open(f'{minim_name}.dat', 'a+') as minimf:
            minimf.write('potential_energy\tconstraint_energy\tconstraint_strength\n')
            for i in range(len(reporter.potential_energies)):
                potential_energy = float(reporter.potential_energies[i])
                constraint_energy = float(reporter.constraint_energies[i])
                constraint_strength = float(reporter.constraint_strengths[i])
                
                minimf.write(f'{potential_energy}\t{constraint_energy}\t{constraint_strength}\n')
