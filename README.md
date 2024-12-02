# md equilibration with OpenMM

This is a script that performs the main equilibration steps of a solvated protein(-ligand) system, prepared with Amber, using OpenMM.
The reason I wrote this is because amber cannot gradually decrease position restraints. Doing multiple runs with discrete restraint weight values is not only cumbersome, it is also not ideal for a smooth system equilibration.