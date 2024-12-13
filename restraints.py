# restraints.py
from openmm import CustomExternalForce
from openmm.unit import Quantity
from simulation import Equilibration
import numpy as np

class PositionRestraints:

    def __init__(self, equilibration:Equilibration, config):
        self.system = self.equilibration.system
        self.equilibration = equilibration
        self.topology = equilibration.prmtop.topology
        self.config = config
        self.equilibration.position_restraints.append(self)

    def apply(self, minim:bool=False):
        
        if minim:
            positions = self.equilibration.inpcrd.positions
            weight = self.config.minim_weight
        else:
            state = self.equilibration.simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
            positions = state.getPositions()
            weight = self.config.weight

        restraint = CustomExternalForce(f'{self.self.config.name}*periodicdistance(x, y, z, x0, y0, z0)^2')
        restraint.addGlobalParameter(self.self.config.name, weight)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

        for atom, position in zip(self.topology.atoms(), positions):
            if self.config.filter_func(atom):
                x0, y0, z0 = position
                restraint.addParticle(atom.index, [x0, y0, z0])

        self.system.addForce(restraint)
        self.restraint = restraint
        self.restraint.setName(self.self.config.name)

    def remove(self):

        for index in range(self.equilibration.system.getNumForces()):
            force = self.equilibration.system.getForce(index)
            if force.getName() == self.restraint.getName():
                self.system.removeForce(index)
                self.equilibration.simulation.context.reinitialize(preserveState=True)


    def update_weight(self, weight:float):

        weight = weight * self.config.weight.unit
        self.equilibration.simulation.context.setParameter(self.config.name, weight)
    

    def get_weight_list(self):

        weight_val = self.config.get_value('weight')
        start_time = self.config.start_time
        stop_time = self.config.stop_time
        dt = self.equilibration.config.dt
        npt_restr_time = self.equilibration.config.npt_restr_time
        
        stage1 = [ weight_val ] * int(np.round(start_time / dt))

        stage2_numsteps = int(np.round( stop_time/dt - start_time/dt ))
        stage2_x = list(range(stage2_numsteps))

        line_params = [weight_val, stage2_numsteps]
        if self.config.decay_func is self.config.exp_decay:
            line_params.append(self.config.decay_rate)
        
        stage2 = [self.config.decay_func(*line_params, x) for x in stage2_x]

        stage3 = [0.0] * int( np.round(npt_restr_time/dt - stop_time/dt) )
    
        with open(f'{self.config.name}_weights.csv', 'a+') as rstfile:
            rstfile.write('weight,name\n')
            rstfile.writelines(f'{str(i)},{self.config.name}\n' for i in [*stage1, *stage2, *stage3])
        
        weight_list = [*stage1, *stage2, *stage3]
        self.weight_list = weight_list