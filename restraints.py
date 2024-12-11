# restraints.py
from openmm import CustomExternalForce
from openmm.unit import Quantity
from simulation import Equilibration

class PositionRestraints:

    def __init__(self, equilibration:Equilibration, posres_name:str, weight:Quantity, filter_func:function):
        self.system = self.equilibration.system
        self.equilibration = equilibration
        self.topology = equilibration.prmtop.topology
        self.positions = equilibration.simulation.state.getPositions()

        self.posres_name = posres_name
        self.filter_func = filter_func

    def apply(self, weight:Quantity) -> int:
        
        restraint = CustomExternalForce(f'{self.posres_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
        restraint.addGlobalParameter(self.posres_name, weight)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

        for atom, position in zip(self.topology.atoms(), self.positions):
            if self.filter_func(atom):
                x0, y0, z0 = position
                restraint.addParticle(atom.index, [x0, y0, z0])

        force_index = self.system.addForce(restraint)
        self.force_index = force_index

    def remove(self):
        self.equilibration.system.removeForce(self.force_index)
        self.equilibration.simulation.context.reinitialize(preserveState=True)

    def update_weight(self, weight:Quantity):
        self.equilibration.simulation.context.setParameter(self.posres_name, weight)