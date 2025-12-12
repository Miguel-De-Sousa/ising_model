PHY2027 Final Project: 2D Ising Model Monte Carlo Simulation
Author: Miguel de Sousa Date: December 2025

📝 Overview
This repository contains the source code and documentation for the PHY2027 Scientific Programming in C final project: a simulation of the Two-Dimensional Ising Model using the Metropolis Monte Carlo algorithm.

The Ising model is a fundamental tool in statistical mechanics used to study phase transitions and critical phenomena in simple magnetic systems.

✨ Key Features & Technical Details
Metropolis Algorithm: Simulates the system's evolution by calculating the change in energy (ΔH) for a spin flip, rather than recalculating the full system energy. This is a crucial optimization for high-speed Monte Carlo simulations.

Memory Optimization: The spin states are stored as char (1 byte) instead of int (4 bytes) to significantly reduce memory footprint for large lattices (N×N).

Periodic Boundary Conditions (PBC): Correctly implemented using modulo arithmetic (% N) to simulate an infinite lattice and ensure accurate nearest-neighbour interactions.

Data Structure: Uses a struct (IsingLattice) to manage all necessary parameters (N, T, B, spin array pointer, energy, magnetisation) for clean and efficient data handling.

Scientific Output: Data is written to a CSV file, ensuring compatibility with standard scientific plotting tools (e.g., Python, Gnuplot).

💡Thank you
Thanks for taking a look at my project!!
