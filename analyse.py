import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md

plt.style.use('seaborn-v0_8')


class Analyse:

    def __init__(self, equilibration):
        self.equilibration = equilibration

    def frame2time(self, frame_list, dec=2):
        write_step = self.equilibration.config.report_interval
        time_list = frame_list * ( 0.002/1000 * write_step )
        return(round(time_list, dec))


    def plot_minimizations(self):

        for i in self.equilibration.minimizations:

            data = pd.read_csv(i+'.dat', sep='\t')
            data = data.astype(float)
            
            plt.figure(tight_layout=True, figsize=[4.0,6.0], dpi=300)
            plt.subplot(211)
            plt.plot(data.index, data['potential_energy'], '-')
            plt.xlabel('Step')
            plt.ylabel('Potential energy (kJ/mol)')
            plt.subplot(212)
            plt.plot(data.index, data['constraint_energy'], '-')
            plt.xlabel('Step')
            plt.ylabel('Constraint energy (kJ/mol)')
            
            plt.savefig(i+'.png')
    

    def calc_rmsds(self, lig_resname=None):

        top = self.equilibration.config.top_name
        last_minim_name = self.equilibration.minimizations[-1]
        trajectories = [i+'.dcd' for i in self.equilibration.simulations]

        traj = md.load(trajectories, top=top)
        ref = md.load(f"{last_minim_name}.rst7", top=top)

        bb_ix = traj.topology.select("protein and backbone")
        ptn_ix = traj.topology.select("protein")

        traj_aln = traj.superpose(ref[0], atom_indices=bb_ix)

        data_dict = {
            'bb_rmsd'  : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=bb_ix), units_in='nanometers', units_out='angstroms'),
            'ptn_rmsd' : md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=ptn_ix), units_in='nanometers', units_out='angstroms'),
            'time'     : self.frame2time(pd.Series(range(1, traj.n_frames+1)))
        }
        if lig_resname is not None:
            lig_ix = traj.topology.select(f"resname {lig_resname}")
            data_dict['lig_rmsd'] = md.utils.in_units_of(md.rmsd(target=traj_aln, reference=ref, atom_indices=lig_ix), units_in='nanometers', units_out='angstroms')

        trajectories = []
        times = [getattr(self.equilibration.config, i) for i in ['heat_time', 'nvt_time', 'npt_restr_time', 'npt_time']]
        time_sums = np.cumsum(times)
        self.time_sums = time_sums
        unit = time_sums[0].unit

        for time in data_dict['time']:
            if time * unit <= time_sums[0]:
                trajectories.append('heat')
            elif time * unit <= time_sums[1]:
                trajectories.append('nvt')
            elif time * unit <= time_sums[2]:
                trajectories.append('npt_restr')
            elif time * unit <= time_sums[3]:
                trajectories.append('npt')

        data = pd.DataFrame.from_dict(data_dict)
        data['traj'] = trajectories
        data.to_csv('equil_rmsd.csv')

        return(data)


    def plot_rmsd(self, lig_resname=None):

        data = pd.read_csv('equil_rmsd.csv', index_col=0)
        unit = self.time_sums[0].unit
        time_sums = [i.value_int_unit(unit) for i in self.time_sums]

        def make_fig():
            plt.figure(tight_layout=True, figsize=[5.0,3.0], dpi=300)
        def make_hlines():
            plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=time_sums[0])
            plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=time_sums[1])
            plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=time_sums[2])
            plt.axvline(linestyle='--', color='#c29ebe', linewidth=1.0,  x=time_sums[3])
        def make_deco():
            plt.legend()
            plt.xlabel('Time (ns)')
            plt.ylabel('RMSD ($\mathrm{\AA}$)')

        make_fig()
        make_hlines()
        plt.plot(data['time'], data['bb_rmsd'], '-', label='backbone', color='#dd842e')
        plt.plot(data['time'], data['ptn_rmsd'], '-', label='not H', color='#519ca2')
        make_deco()
        plt.savefig('rmsd.png')

        if lig_resname is not None:

            make_fig()
            make_hlines()
            plt.plot(data['time'], data['lig_rmsd'], '-', label='not H', color='#519ca2')
            make_deco()
            plt.savefig('lig_rmsd.png')


    # def plot_restr_weights(self):

    #     weights = {}
    #     for i in self.equilibration.position_restraints:
    #         name = i.config.name
    #         weights.append(name, i.weight_list)

    #     filetypes = {
    #         'ptn_restr':[file for file in files if "lig" not in str(file.stem)],
    #         'lig_restr':[file for file in files if "lig" in str(file.stem)]
    #     }
        
    #     for filetype in filetypes:

    #         fig = plt.figure(tight_layout=True, figsize=[5.0,2.0], dpi=300)
    #         colors = ['#dd842e','#519ca2']
    #         color_ix = 0
    #         for file in filetypes[filetype]:
    #             data = pd.read_csv(file)
    #             data['index'] = data.reset_index()['index']+1
    #             data['time'] = frame2time(data['index'], dec=10, write_step=1)
    #             plt.plot(data['time'], data['weight'], label=file.stem, color=colors[color_ix])
    #             color_ix+=1
            
    #         plt.xlabel('Time (ns)')
    #         plt.ylabel('RWT (kcal/mol/$\mathrm{\AA}^2$)')
    #         plt.legend()

    #         plt.savefig(filetype + '.png')