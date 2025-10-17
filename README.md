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

For an MD project, it is good practice to have note of the protocol and all the parameters used. Since this program encodes the whole protocol and MD parameters used in an equilibration run, I developed it not as a standalone application, but as a folder that resides in your simulation project. You can therefore clone the repo and edit the parameters for your use case in that specific project, thus making it future-proof: you will always be able to see exactly what the protocol and parameters used were, even if the code updates on GitHub or you transfer the project folder between computers.
Example project folder:

```
md_simulation/
├── equilibration
│   ├── MDProtocol
│   ├── rep_01
│   ├── rep_02
│   └── rep_03
├── protein.pdb
├── system.prmtop
├── system.rst7
└── tleap.in
```

Note that, because this was developed as a python package, the `MDProtocol` folder must be on the current working directory. But for replicates, since the configuration between them is the same, there is no need to copy the `MDProtocol` folder inside each one, simply create a symbolic link in each replicate folder to the `MDProtocol` one folder up. Thus:

```
├── equilibration
│   ├── MDProtocol
│   ├── rep_01
│   │   └── MDProtocol -> ../MDProtocol
```

This program requires the python packages OpenMM, Parmed, and MDTraj, ideally in a clean virtual environment.

Then, simply run with:

```
python -m MDProtocol
```

### Configuration

(Work in progress)

All the configuration options the large majority of users will ever need are at the file `__main__.py`. We will break down this file, step by step. 

First the system and simulation options are configured. This block of code defines the simulation parameters:
```python
simconf = SimulationConfig(
    # add here all options that should not be kept default
    top_name='complex.prmtop',
    crd_name='complex.rst7',
)
```

The default options are as follows (and can also be found in `config.py`):
```py
top_name:         str   = 'system.prmtop'
crd_name:         str   = 'system.rst7'

dt:               float = 0.002
tau_t:            float = 1.0
cutoff:           float = 10.0
report_interval:  int   = 5000

heat_time:        float = 3.0   # total time for heating nvt
nvt_time:         float = 2.0   # total time for nvt at target temp
npt_restr_time:   float = 15.0  # total time for npt with gradually reduced restraints
npt_time:         float = 5.0   # total time for npt with no restraints

start_temp:       float = 100.0
end_temp:         float = 300.0

pressure:         float = 1.0
```

Thus, if the user wishes to change any of these options, they just need to alter `simconf`, for example:
```python
simconf = SimulationConfig(
    # add here all options that should not be kept default
    top_name='complex.prmtop',
    crd_name='complex.rst7',
    report_interval=2500
)
```

Next, ligand information, structural waters and position restraints are defined. If you have a ligand in your system, provide its residue names here:
```py
lig_resname = ['LIG'] # IMPORTANT: if no ligand, use lig_resname = []
```
If your ligand contains multiple parts with different residue names, add all of them to the list, example `['PT1','PT2']`
If you do not have a ligand, keep:
```py
lig_resname = [] # IMPORTANT: if no ligand, use lig_resname = []
```
Then, provide the atom names of the anchor portion of the ligand (the program will automatically identify the extension portion):
```py
lig_anchor_atoms = ['C1', 'C2', 'C3', 'C4', 'C5', 'O3', 'C9', 'C12', 'C13', 'C14', 'C15', 'O10', 'C20', 'C21', 'C22', 'C23', 'C24', 'O16', 'C28', 'C31', 'C32', 'C33', 'C34', 'O22']
```
If you have structural waters, provide their residue numbers as strings (within quotation marks):
```py
structural_waters = ['163', '164', '165', '166'] # residue numbers
```
If you don't have any, comment this line out or delete it.

Finally, we can configure the position restraints. Each posres block defines these configurations:
- `name`: the posres' name. Default "posres_bb"
- `weight`: its weight in kilocalories/mole/angstroms**2, default 10.0.
- `minim_weight`: its weight in the minimization stage (this way one weight can be used for minimizing and another for equilibrating). Default 10.0.
- `start_time`: the equilibration time point where its weight starts to be decreased (ns). Default 5.0.
- `end_time`: the equilibration time point where its weight reaches 0. Default 15.0.
- `lig_resname`: the residue names of the ligand as a list of strings. Default is None.
- `lig_anchor_atoms`: the ligand anchor atoms. Default None.
- `structural_waters`: the list of structural waters. Default None.
- `mask_func_name`: the name of the mask function that this specific position restraint will use. It can be: "posres_bb_mask", "posres_sc_mask", "posres_liganc_mask", "posres_ligext_mask", or "posres_water_mask". Default "posres_bb_mask".
- `decay_func`: use linear decay or exponential decay. Uses exponential by default.
- `decay_rate`: the rate of decay. Default 5.0.

The parameters that are to be kept at default value do not need to be explicitly defined.
So let's start defining the config for the backbone position restraints. We will put it in the variable `config_posres_bb`. The restraint's name will be `"posres_bb"`. And since this is the backbone postres, the atom mask it will use is the `"posres_bb_mask"`. We will need to add our ligand specification. We will keep all the other values as default and therefore omit them. Thus:

```py
config_posres_bb = RestraintConfig( # protein backbone
    name      = 'posres_bb',
    mask_func_name = 'posres_bb_mask',
    lig_resname=lig_resname
)
```

Next, we configure the posres for the protein sidechain. It will be called `"posres_sc"`:

```py
config_posres_sc = RestraintConfig( # protein sidechain
    name       = 'posres_sc',
    lig_resname=lig_resname,
    weight     = 5.0, # change the weight value from 10.0 to 5.0
    mask_func_name  = 'posres_sc_mask', # since this is posres for sidechain, the mask used will be the one for sidechain
    start_time = 0, # the sidechain will start being release at the start of npt
    end_time   = simconf.get_value('npt_restr_time') / 3 # it will reach 0 at total npt time divided by 3. Similarly, you could just say `end_time = 5.0` if your npt time is 15.0
)
```

For configuring the posres on the ligand anchor, according to the protocol the options should follow the ones set for the protein backbone, but we should give it the information we set previously in the variables `lig_resname` and `lig_anchor_atoms`, so:

```py
config_posres_liganch = RestraintConfig( # ligand anchor
    name        = 'posres_liganc',
    mask_func_name   = 'posres_liganc_mask',
    lig_resname=lig_resname,
    lig_anchor_atoms=lig_anchor_atoms
)
```

Similarly the other restraints are configured:

```py
config_posres_ligext = RestraintConfig( # ligand extension
    name = 'posres_ligext',
    weight = 5.0, # same weight as the protein sidechain
    mask_func_name = 'posres_ligext_mask', # use the corresponding atom mask
    lig_resname=lig_resname, # provide ligand information set previously
    lig_anchor_atoms=lig_anchor_atoms, # provide ligand information set previously
    start_time = 0, # like the protein sidechain, starts decreasing at the start of the npt simulation 
    end_time = simconf.get_value('npt_restr_time') / 3 # and reaches 0 at a third of the simulation
)
config_posres_wat = RestraintConfig( # structural waters
    name = 'posres_wat',
    weight = 1.0, 
    mask_func_name = 'posres_water_mask',
    structural_waters = structural_waters, # give the structural waters' resnumbers defined previously
    start_time = 0,
    end_time = simconf.get_value('npt_restr_time') / 3
)
```

You can comment out or delete any of these blocks if you will not use the corresponding posres.
The next step only initialises the position restraint objects within the program using the configurations we just defined, but they do not get applied to the simualation itself yet:

```py
# start system
equilibration = Equilibration(simconf)

# initialize position restraints (this does not apply them yet)
posres_bb = PositionRestraints(equilibration, config_posres_bb)
posres_sc = PositionRestraints(equilibration, config_posres_sc)
posres_liganc = PositionRestraints(equilibration, config_posres_liganch)
posres_ligext = PositionRestraints(equilibration, config_posres_ligext)
posres_waters = PositionRestraints(equilibration, config_posres_wat)
# ^ comment out or delete any of the lines corresponding to posres that will not be used
```

Now we can apply the position restraints that will be used for minimization (delete any that will not be used in your protocol):
```py
# apply all restraints except struct waters
posres_bb.apply(minim=True)
posres_sc.apply(minim=True)
posres_liganc.apply(minim=True)
posres_ligext.apply(minim=True)
```
Note that we are applying the minimization weights, thus we used `minim=True`
Then we can minimize the system:
```py
equilibration.minimize('minim1')
```

Next, we will remove the ligand extension and protein sidechain posres so that the system can minimize again only keeping backbone and anchor restrained:

```py
# remove posres sc and ligext:
posres_sc.remove()
posres_ligext.remove()
# minimize 2nd round
equilibration.minimize('minim2')
```

And then we remove the remaining restraints to perform a free minimization:

```py
# remove posres bb and liganc:
posres_bb.remove()
posres_liganc.remove()
# minimize 3rd round:
equilibration.minimize('minim3')
```

Next, the nvt simulation can be started. We apply SHAKE and then reapply the position restraints, now with the simulation weights (not the minimization weights):

```py
# start SHAKE
equilibration.apply_hbond_constraints()

# apply all position restraints for simulation
for restr in equilibration.position_restraints:
    restr.apply()
```

Next we run nvt while increase the temperature, followed by nvt at target temperature:
```py
# increase temp linearly
equilibration.heat('nvt_heat')
# keep nvt at target temperature
equilibration.nvt('nvt')
```

We add the barostat and then run the npt stage where all the restraints are gradually reduced to 0 obeying the settings and timings we defined earlier:

```py
add_barostat(simconf, equilibration)
equilibration.npt_posres('npt_posres')
```

Then we run a quick npt without any position restraint applied:

```py
# even though restraint weights got to 0, I remove them
# to avoid unnecessery calculations:
for restr in equilibration.position_restraints:
    restr.remove()

equilibration.npt('npt')
```

Finally, the program will generate some nice plots automatically:
```py
analysis = Analyse(equilibration=equilibration)

analysis.calc_rmsds(lig_resname=lig_resname)
analysis.plot_minimizations()
analysis.plot_rmsd(lig_resname=lig_resname)
analysis.plot_restr_weights()
```

Delete any line corresponding to an analysis you do not want done.


### Acknowledgements

Thanks [Max](https://orcid.org/0009-0008-7252-9348), for the great discussions about this protocol.
