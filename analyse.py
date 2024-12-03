import matplotlib.pyplot as plt
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
    
plot_minimizations()