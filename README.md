# MD equilibration with OpenMM

This is a script that performs the main equilibration steps of a solvated protein(-ligand) system, prepared with Amber, using OpenMM.
The reason I wrote this is because amber cannot gradually decrease position restraints. Doing multiple runs with discrete restraint weight values is not only cumbersome, it is also not ideal for a smooth system equilibration.

# Protocol

PME, with 10A cutoff.
Position restraints:
- backbone: 10kcal/mol/A**2
- sidechain, not hydrogens: 5kcal/mol/A**2

## Minimization

- Minim1: minimize with both restraints
- Minim2: minimize only with bb restraint
- Minim3: minimize with no restraints

## NVT

Using NoseHoover thermostat.

- Slowly heat the system from 100K to 300K, for 3ns
- Keep system at 300K for 2ns

## NPT

Turn on MonteCarlo barostat.

- For 15ns, slowly release each restraint:
    - during the first third, gradually reduce sidechain restraints
    - during the last two thirds, gradually reduce brackbone restraints
- For 5ns, simulate NPT with no restraints