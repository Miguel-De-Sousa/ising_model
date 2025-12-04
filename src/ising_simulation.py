import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

df = pd.read_csv('src/ising_data.csv', header=None)

steps = df.iloc[:, 0].values
energies = df.iloc[:, 1].values
mags = df.iloc[:, 2].values
spin_data = df.iloc[:, 3:].values         


grid_size = 512
lattices = spin_data.reshape(-1, grid_size, grid_size) 

fig, ax = plt.subplots(figsize=(8, 8)) 
im = ax.imshow(lattices[0], cmap='gray', vmin=0, vmax=1, interpolation='nearest')
ax.set_title(f'Step {steps[0]:6d} | E = {energies[0]:7.1f} | M = {mags[0]:+.3f}',
             fontsize=16, pad=20, family='monospace')
ax.axis('off')

def update(frame):
    im.set_array(lattices[frame])
    ax.set_title(f'Step {steps[frame]:6d} | E = {energies[frame]:7.1f} | M = {mags[frame]:+.3f}',
                 fontsize=16, pad=20, family='monospace')
    return [im]

ani = FuncAnimation(fig, update, frames=len(steps), interval=100, blit=True, repeat=False)
ani.save('src/ising_evolution.gif',
         writer=PillowWriter(fps=120),
         dpi=100)                      

plt.close(fig)