# restraints.py
from openmm import CustomExternalForce, System
from openmm.unit import kilocalories_per_mole, angstroms, Quantity
from typing import Tuple

class PositionRestraints:
    def __init__(self, system:System, topology, positions, config):
        self.system = system
        self.topology = topology
        self.positions = positions
        self.config = config
        self.restraints = {}

    def apply_restraints(self, restraint_type: str, weight:Quantity, ligand_center_atoms=None) -> int:
        """
        Applies a position restraint based on the type.
        """
        
        if restraint_type == 'bb':
            filter_func = self._filter_bb
            posres_name = 'k_bb'
        elif restraint_type == 'sc':
            filter_func = self._filter_sc
            posres_name = 'k_sc'
        elif restraint_type == 'lig':
            filter_func = lambda atom: self._filter_lig(atom, ligand_center_atoms)
            posres_name = 'k_lig'
        elif restraint_type == 'lig_center':
            filter_func = lambda atom: self._filter_lig_center(atom, ligand_center_atoms)
            posres_name = 'k_ligc'
        else:
            raise ValueError(f"Unknown restraint type: {restraint_type}")

        restraint = CustomExternalForce(f'{posres_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
        restraint.addGlobalParameter(posres_name, weight)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

        for atom, position in zip(self.topology.atoms(), self.positions):
            if filter_func(atom):
                x0, y0, z0 = position
                restraint.addParticle(atom.index, [x0, y0, z0])

        force_index = self.system.addForce(restraint)
        self.restraints[restraint_type] = force_index
        return force_index

    def _filter_bb(self, atom):
        return (
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and
            atom.name in ['C', 'CA', 'N']
        )

    def _filter_sc(self, atom):
        return (
            atom.residue.name not in ['HOH', 'WAT', 'Na+', 'Cl-', 'LIG'] and
            atom.element.symbol != "H" and
            atom.name not in ['C', 'CA', 'N']
        )

    def _filter_lig(self, atom, ligand_center_atoms):
        return (
            atom.residue.name == 'LIG' and
            atom.element.symbol != "H" and
            atom.name not in ligand_center_atoms
        )

    def _filter_lig_center(self, atom, ligand_center_atoms):
        return (
            atom.residue.name == 'LIG' and
            atom.name in ligand_center_atoms
        )
