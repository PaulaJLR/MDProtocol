## Automated protein(-ligand) equilibration with OpenMM to overcome the pitfalls of common protocols with a reliable, modern approach.

In MD simulation protocols, proper equilibration is an essential step to ensure the system exhibits biologically representative dynamics, and thus to prepare the system for accurate data collection.

This repository implements equilibration steps of a solvated protein system, prepared with Amber, using OpenMM. The protocol implemented here was designed to approach equilibration in the smoothest way possible. To this end, it addresses significant position restraint limitations that are currently not easy to overcome in common MD software (e.g. Amber, Gromacs, or NAMD).

### The position restraint problem

The first step to an MD equilibration is to gradually heat the system using the NVT ensemble. Then, once at target temperature, we turn on the basrostat for NPT equilibration. This means that the solute will be subjected to non physiological conditions and forces, which could lead to unrepresentative conformational changes. To mitigate this, it is common to implement cartesian position restraints (posres) on the solute, so that it is kept near the starting conformation and is not allowed to change much during this initial stage.

However, few software solutions allow gradual removal of the posres. As a result, protocols adopt a "staircase" approach, in which the user is responsible for scripting multiple MD runs where the restraint weight is reduced to 0 in descrete steps. In this case, it is extremely cumbersome to implement numerous steps. Necessarily, the steps are far enough from each other that each decrease in restraint weight will introduce a very sharp change to the system. This too is undesirable, since it can also lead to unrepresentative conformational changes.

Another problem in the staircase solution is that it is common practice to treat backbone and sidechains separately, but it is not always straightforward to define two separate atom groups. 

### This approach

I developed this code to implement gradual position restraint reduction to my Amber equilibration protocol, using OpenMM. 
Its main advantage is that it implements **exponential decay** of position restraint weights and, importantly, the **weights are updated at every single simulation step** (without significant performance loss), resulting in extremely smooth restraint reduction. In addition, side chain and backbone groups are treated separately, such that sidechain posres can be reduced to 0 (or any other weight), before the backbone posres can start decreasing.

Initally, the idea was to only implement linear restraint weight reduction (which is still supported by this script). However, it was apparent that, in a relatively stable protein system, the lower restraint values accounted for most of the RMSD changes, while the higher values still maintained a very narrow RMSD on the low range. Therefore, an exponential decay on the weight values ensures that we reduce the higher values more quickly, and the lower values even more gradually. This results in a very satisfying, smoothly increasing RMSD plot.

(example images soon)

### The Protocol

Out of the box, his script support systems with any combination of these components: protein, ligand, structural waters.
For applying position restraints, the protein is divided into backbone and sidechain groups. Similarly, the ligand may be divided into anchor (central and/or more rigid scaffold) and extension.
Most of the protocol parameters, such as the inclusion of posres and their weights, minimization rounds, system temperature, etc. are easily configurable; please see "usage" section.

PME, with 10A cutoff.

Position restraints for initial minimization (kcal/mol/A**2):
- backbone: 10
- sidechain and not hydrogens: 10
- ligand anchor: 10
- ligand extension: 10
- structural waters: none

Position restraints for simulation (kcal/mol/A**2):
- backbone: 10
- sidechain and not hydrogens: 5
- ligand anchor: 10
- ligand extension: 5
- structural waters: 1


#### Minimization

It is by default done in three rounds:

- 1. minimize with all restraints except structural waters
- 2. minimize only with protein backbone and ligand anchor restraints
- 3. minimize with no restraints

#### NVT

Parameters:
- dt = 2fs
- NoseHoover thermostat with tau_t = 1ps
- apply HBond length constraints

Simulation:
- Slowly heat the system from 100K to 300K, for 3ns
- Keep system at 300K for 2ns

#### NPT

Turn on MonteCarlo barostat.

- Gradually reduce position restraints over a 15ns period using an exponential decay schedule:
    - first third (0-5ns), gradually reduce protein sidechain and ligand extension restraints
    - last two thirds (5-15ns), gradually reduce protein backbone and ligand anchor restraints
- For additional 5ns, perform unrestrained NPT simulation.

#### Analyses

This script automatically performs RMSD analyses over all equilibration stages and generates the following plots:

- Minimization energies
- RMSD evolution for protein (backbone and sidechain)
- RMSD evolution for ligand
- Restraint weights applied throughout the NPT with posres reduction

### Usage

(Coming soon)

### Acknowledgements

Thanks [Max](https://orcid.org/0009-0008-7252-9348), for the great discussions about this protocol.