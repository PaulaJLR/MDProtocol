import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
from equilibrate import heat_time, nvt_time, npt_restr_time, npt_time
from openmm.unit import *

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
    

def frame2time(frame_list, dec=2, write_step=5000):
    time_list = frame_list * ( 0.002/1000 * write_step )
    return(round(time_list, dec))


def calc_rmsds():

    

    top = "complex.prmtop"
    trajectories = ["nvt_heat.dcd", "nvt.dcd", "npt_posres.dcd", "npt1.dcd"]

    traj = md.load(trajectories, top=top)
    ref = md.load("minim3.rst7", top=top)

    bb_ix = traj.topology.select("protein and backbone")
    ptn_ix = traj.topology.select("protein")
    lig_ix = traj.topology.select("resname LIG")

    traj_aln = traj.superpose(ref[0], atom_indices=bb_ix)

    data_dict = {
        'bb_rmsd'  : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=bb_ix), units_in='nanometers', units_out='angstroms'),
        'ptn_rmsd' : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=ptn_ix), units_in='nanometers', units_out='angstroms'),
        'lig_rmsd' : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=lig_ix), units_in='nanometers', units_out='angstroms'),
        'time'     : frame2time(pd.Series(range(1, traj.n_frames+1)))
    }

    trajectories = []
    for time in data_dict['time']:
        if time * nanoseconds <= heat_time:
            trajectories.append('heat')
        elif time * nanoseconds <= nvt_time + heat_time:
            trajectories.append('nvt')
        elif time * nanoseconds <= npt_restr_time + heat_time + nvt_time:
            trajectories.append('npt_restr')
        elif time * nanoseconds <= npt_time + npt_restr_time + heat_time + nvt_time:
            trajectories.append('npt')

    data = pd.DataFrame.from_dict(data_dict)
    data['traj'] = trajectories
    data.to_csv('equil_rmsd.csv')

    return(data)


def plot_rmsds():

    data = pd.read_csv('equil_rmsd.csv', index_col=0)

    fig = plt.figure(tight_layout=True, figsize=[5.0,3.0], dpi=300)

    plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=heat_time.value_in_unit(nanoseconds))
    plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=(nvt_time + heat_time).value_in_unit(nanoseconds))
    plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=(npt_restr_time + nvt_time + heat_time).value_in_unit(nanoseconds))
    plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=(npt_time + npt_restr_time + nvt_time + heat_time).value_in_unit(nanoseconds))

    plt.plot(data['time'], data['bb_rmsd'], '-', label='backbone', color='#dd842e')
    plt.plot(data['time'], data['ptn_rmsd'], '-', label='not H', color='#519ca2')

    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD ($\mathrm{\AA}$)')
    plt.savefig('rmsd.png')


# plot_minimizations()
# calc_rmsds()
plot_rmsds()
