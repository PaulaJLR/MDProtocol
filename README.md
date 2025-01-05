## Automated protein(-ligand) equilibration with OpenMM to overcome the pitfalls of common protocols with a reliable, modern approach.

In MD simulation protocols, proper equilibration is an essential step to ensure the system exhibits biologically representative dynamics, and thus to prepare the system for accurate data collection.

This repository implements equilibration steps of a solvated protein system, prepared with Amber, using OpenMM. The protocol implemented here was designed to approach equilibration in the smoothest way possible. To this end, it addresses significant position restraint limitations that are currently not easy to overcome in common MD software (e.g. Amber, Gromacs, or NAMD).

### The position restraint problem

The first step to an MD equilibration is to gradually heat the system using the NVT ensemble. Then, once at target temperature, we turn on the basrostat for NPT equilibration. This means that the solute will be subjected to non physiological conditions and force changes, which could lead to unrepresentative conformational changes. To mitigate this, it is common to implement cartesian position restraints on the solute, so that it is kept near the starting conformation and is not allowed to change much during this initial stage.

However, few software solutions allow gradual removal of the position restraints. As a result, protocols adapt a "staircase" approach, in which the user is responsible for scripting multiple MD runs where the restraint weight is reduced to 0 in descrete steps. In this case, it is extremely cumbersome to implement numerous steps. Necessarily, the steps are far enough from each other that each decrease in restraint weight will introduce a very sharp change to the system. This is also undesirable, since it can also lead to unrepresentative conformational changes.

### This approach

I developed this code to implement gradual position restraint reduction to my Amber MD equilibration protocol. Its main advantage is that it implements **exponential decay** of position restraint weights and, importantly, the **weights are updated at every single simulation step**, resulting in extremely smooth reduction.

Initally, the idea was to only implement linear restraint weight reduction (which is still supported by this script). However, it was apparent that, in a relatively stable protein system, the lower restraint values accounted for most of the RMSD changes, while the higher values still maintained a very narrow RMSD on the low range. Therefore, an exponential decay on the weight values ensures that we reduce the higher values more quickly, and the lower values even more gradually. This results in a very satisfying, smoothly increasing RMSD plot.

### The Protocol

PME, with 10A cutoff.
Position restraints:
- backbone: 10kcal/mol/A**2
- sidechain, not hydrogens: 5kcal/mol/A**2

#### Minimization

- Minim1: minimize with both restraints
- Minim2: minimize only with bb restraint
- Minim3: minimize with no restraints

#### NVT

Using NoseHoover thermostat.

- Slowly heat the system from 100K to 300K, for 3ns
- Keep system at 300K for 2ns

#### NPT

Turn on MonteCarlo barostat.

- For 15ns, slowly release each restraint:
    - during the first third, gradually reduce sidechain restraints
    - during the last two thirds, gradually reduce brackbone restraints
- For 5ns, simulate NPT with no restraints

### Acknowledgements

Thanks [Max](https://orcid.org/0009-0008-7252-9348), for the great discussions about this protocol.