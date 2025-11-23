import numpy as np 
import matplotlib.pyplot as plt

def plot_ising_data(file_path):
    data = np.loadtxt(file_path, delimiter=',', skiprows=1)
    steps = data[:, 0]
    energies = data[:, 1]
    magnetisations = data[:, 2]

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Steps')
    ax1.set_ylabel('Energy', color=color)
    ax1.plot(steps, energies, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  
    color = 'tab:blue'
    ax2.set_ylabel('Magnetisation', color=color)  
    ax2.plot(steps, magnetisations, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  
    plt.title('2D Ising Model Simulation Data')
    plt.show()

if __name__ == "__main__":
    plot_ising_data('src/ising_data.csv')
    