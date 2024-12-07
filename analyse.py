import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md

plt.style.use('seaborn-v0_8')

def plot_minimizations():

    for i in ['minim1', 'minim2', 'minim3']:

        energies = open(i+'.dat', 'r').readlines()
        energies = [float(i.strip()) for i in energies]
        
        fig = plt.figure(tight_layout=True, figsize=[4.0,3.0], dpi=300)
        plt.plot(list(range(len(energies))), energies, '-')
        plt.xlabel('Step')
        plt.ylabel('System energy (kJ/mol)')
        plt.savefig(i+'.png')
    

def calc_rmsds():

    top = "complex.prmtop"
    trajectories = ["nvt_heat.dcd", "nvt.dcd", "npt_posres.dcd", "npt1.dcd"]

    traj = md.load(trajectories, top=top)
    ref = md.load("minim3.rst7", top=top)

    bb_ix = traj.topology.select("protein and backbone")
    lig_ix = traj.topology.select("resname LIG")

    traj_aln = traj.superpose(ref[0], atom_indices=bb_ix)

    data_dict = {
        'protein_rmsd' : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=bb_ix), units_in='nanometers', units_out='angstroms'),
        'lig_rmsd'     : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=lig_ix), units_in='nanometers', units_out='angstroms')
    }
    data = pd.DataFrame.from_dict(data_dict)
    data.to_csv('equil_rmsd.csv')

    return(data)

# plot_minimizations()
calc_rmsds()