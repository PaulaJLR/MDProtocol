import numpy as np
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds
from openmm import NoseHooverIntegrator
from openmm.unit import kelvin
import parmed as pmd
from .reporters import MyMinimizationReporter, add_reporters
from .config import SimulationConfig
from .tools import save_rst7

class Equilibration:

    def __init__(self, config:SimulationConfig):
        self.config = config
        
        self.prmtop = AmberPrmtopFile(config.top_name)
        self.inpcrd = AmberInpcrdFile(config.crd_name)
        self.system = self.prmtop.createSystem(
            nonbondedMethod=PME,
            nonbondedCutoff=self.config.cutoff,
            constraints=None
        )
        self.integrator = NoseHooverIntegrator(
            self.config.start_temp,
            self.config.tau_t,
            self.config.dt
        )
        self.simulation = Simulation(
            self.prmtop.topology,
            self.system,
            self.integrator
        )
        self.simulation.context.setPositions(self.inpcrd.positions)

        self.position_restraints = []
        self.minimizations = []
        self.simulations = []


    def minimize(self, minim_name):

        reporter = MyMinimizationReporter()
        self.simulation.minimizeEnergy(reporter=reporter)
        save_rst7(self, f'{minim_name}.rst7')

        self.minimizations.append(minim_name)
        with open(f'{minim_name}.dat', 'a+') as minimf:
            minimf.write('potential_energy\tconstraint_energy\tconstraint_strength\n')
            for i in range(len(reporter.potential_energies)):
                potential_energy = float(reporter.potential_energies[i])
                constraint_energy = float(reporter.constraint_energies[i])
                constraint_strength = float(reporter.constraint_strengths[i])
                
                minimf.write(f'{potential_energy}\t{constraint_energy}\t{constraint_strength}\n')

    def apply_hbond_constraints(self):

        prev_state = self.simulation.context.getState(getPositions=True)

        self.system = self.prmtop.createSystem(
            nonbondedMethod=PME,
            nonbondedCutoff=self.config.cutoff,
            constraints=HBonds
        )
        self.integrator = NoseHooverIntegrator(
            self.config.start_temp,
            self.config.tau_t,
            self.config.dt
        )
        self.simulation = Simulation(
            self.prmtop.topology,
            self.system,
            self.integrator
        )
        self.simulation.system = self.system
        self.simulation.context.setPositions(prev_state.getPositions())


    def heat(self, heat_name):

        add_reporters(self, self.config, heat_name)

        nsteps = int(np.round(self.config.heat_time / self.config.dt))
        temps = np.linspace(self.config.start_temp.value_in_unit(kelvin), self.config.end_temp.value_in_unit(kelvin), nsteps)
        
        # run
        for i in range(len(temps)):
            temp = temps[i]
            self.integrator.setTemperature(temp * kelvin)
            self.simulation.step(1)

        save_rst7(self, f'{heat_name}.rst7')
        self.simulations.append(heat_name)
    

    def nvt(self, nvt_name):

        add_reporters(self, self.config, nvt_name)
        nsteps = int(np.round(self.config.nvt_time / self.config.dt))

        # run
        self.simulation.step(nsteps)

        save_rst7(self, f'{nvt_name}.rst7')
        self.simulations.append(nvt_name)
    

    def npt_posres(self, npt_posres_name):

        add_reporters(self, self.config, npt_posres_name)
        nsteps = int(np.round(self.config.npt_restr_time / self.config.dt))

        # get weight lists for each restraint
        for restraint in self.position_restraints:
            restraint.get_weight_list()

        # run
        for i in range(nsteps):

            for restraint in self.position_restraints:

                weight = restraint.weight_list[i]
                restraint.update_weight(weight)

            self.simulation.step(1)

        save_rst7(self, f'{npt_posres_name}.rst7')
        self.simulations.append(npt_posres_name)


    def npt(self, npt_name):

        add_reporters(self, self.config, npt_name)
        nsteps = int(np.round(self.config.npt_time / self.config.dt))

        self.simulation.step(nsteps)

        save_rst7(self, f'{npt_name}.rst7')
        self.simulations.append(npt_name)
