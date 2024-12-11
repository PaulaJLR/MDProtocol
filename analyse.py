import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from pathlib import Path
from equilibrate import heat_time, nvt_time, npt_restr_time, npt_time
from openmm.unit import *

plt.style.use('seaborn-v0_8')

def plot_minimizations():

    for i in ['minim1', 'minim2', 'minim3']:

        data = pd.read_csv(i+'.dat', sep='\t')
        data = data.astype(float)
        
        fig = plt.figure(tight_layout=True, figsize=[4.0,6.0], dpi=300)
        plt.subplot(211)
        plt.plot(data.index, data['potential_energy'], '-')
        plt.xlabel('Step')
        plt.ylabel('Potential energy (kJ/mol)')
        plt.subplot(212)
        plt.plot(data.index, data['constraint_energy'], '-')
        plt.xlabel('Step')
        plt.ylabel('Constraint energy (kJ/mol)')
        
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


def plot_restr_weights():

    files = list(Path.cwd().glob("**/*weights.csv"))
    filetypes = {
        'ptn_restr':[file for file in files if "lig" not in str(file.stem)],
        'lig_restr':[file for file in files if "lig" in str(file.stem)]
    }
    
    for filetype in filetypes:

        fig = plt.figure(tight_layout=True, figsize=[5.0,2.0], dpi=300)
        colors = ['#dd842e','#519ca2']
        color_ix = 0
        for file in filetypes[filetype]:
            data = pd.read_csv(file)
            data['index'] = data.reset_index()['index']+1
            data['time'] = frame2time(data['index'], dec=10, write_step=1)
            plt.plot(data['time'], data['weight'], label=file.stem, color=colors[color_ix])
            color_ix+=1
        
        plt.xlabel('Time (ns)')
        plt.ylabel('RWT (kcal/mol/$\mathrm{\AA}^2$)')
        plt.legend()

        plt.savefig(filetype + '.png')


plot_minimizations()
calc_rmsds()
plot_rmsds()
plot_restr_weights()