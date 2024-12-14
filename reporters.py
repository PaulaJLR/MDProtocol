from openmm import MinimizationReporter
from openmm.app import StateDataReporter, CheckpointReporter, DCDReporter

class MyMinimizationReporter(MinimizationReporter):

    # adapted from https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html

    potential_energies = []
    constraint_energies = []
    constraint_strengths = []

    # you must override the report method and it must have this signature.
    def report(self, iteration, x, grad, args):
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


def add_reporters(equilibration, config, out_name):

    equilibration.simulation.reporters.clear()

    equilibration.simulation.reporters.append(StateDataReporter(
        f"{out_name}.log", config.report_interval, step=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=config.total_steps, separator='\t'
    ))
    equilibration.simulation.reporters.append(DCDReporter(f'{out_name}.dcd', config.report_interval))  # Trajectory output
    equilibration.simulation.reporters.append(CheckpointReporter(f'{out_name}.chk', config.report_interval))  # Checkpointing
