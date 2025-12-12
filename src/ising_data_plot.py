import pandas as pd
import matplotlib.pyplot as plt

def plot_ising_data(file_path):
    try:
        df = pd.read_csv(
            file_path, 
            header=None, 
            usecols=[0, 1, 2],
            names=['Steps', 'Energy', 'Magnetisation']
        )
    except FileNotFoundError:
        print(f"Error: File not found at path: {file_path}")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    

    axes[0].plot(df['Steps'], (df['Energy']/(10**4)), label='Total Energy (E)', color='darkred', linewidth=0.8)
    axes[0].set_title(f'Ising Model Evolution')
    axes[0].set_ylabel(r'Total Energy $\times 10^{4}$')
    axes[0].grid(True, linestyle='--', alpha=0.6)
    axes[0].legend(loc='upper right')
    
    axes[1].plot(df['Steps'], df['Magnetisation'], label=r'Magnetisation ($\mid$M$\mid$)', color='darkblue', linewidth=0.8)
    axes[1].set_xlabel('Monte Carlo Steps', labelpad=10)
    axes[1].set_ylabel('Normalised Magnetisation')
    
    axes[1].axhline(0, color='gray', linestyle='-', linewidth=0.7)
    
    axes[1].set_ylim(-1.1, 1.1) 
    axes[1].grid(True, linestyle='--', alpha=0.6)
    axes[1].legend(loc='upper right')
    
    plt.xlim(0, 100000)
    plt.tight_layout()
    plt.savefig("evolution_plot.png")
    plt.show()

data_file = './data/ising_data.csv' 

plot_ising_data(data_file)