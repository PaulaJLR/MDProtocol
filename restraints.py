# restraints.py
from openmm import CustomExternalForce, System
from openmm.unit import kilocalories_per_mole, angstroms, Quantity
from typing import Tuple
from simulation import Equilibration
from config import RestraintConfig

class PositionRestraints:
    def __init__(self, system:System, topology, positions, restraint_config:RestraintConfig, equilibration:Equilibration, posres_name:str, filter_func:function):
        self.system = system
        self.topology = topology
        self.positions = positions
        self.restraint_config = restraint_config
        self.equilibration = equilibration

        self.posres_name = posres_name
        self.filter_func = filter_func

        self.equilibration.position_restraints.append(self)

    def apply_restraints(self, weight:Quantity) -> int:
        """
        Applies a position restraint based on the type.
        """
        
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

        
    def update_weight(self, weight:Quantity):

        self.equilibration.simulation.context.setParameter(self.posres_name, weight)