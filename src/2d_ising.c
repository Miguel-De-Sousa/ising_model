/*
 * Assessment 3 for PHY2027
 * Author: Miguel de Sousa
 * Date: 12/12/25
*/

/*
 * A program to create a general 2-Dimensional Ising model simulation using snapshot frames, converted into video using Python.
 * The program allow the user to select pre-determined or custom parameters for simulation run.
 * The model uses a random generation of up, down spin values in an N x N lattice.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef char Spin; // Define Spin type as char for memory efficiency
static const Spin SPIN_UP   =  1;
static const Spin SPIN_DOWN = -1;

typedef struct { // Structure to hold the Ising lattice data
    int N;
    Spin *spin;
    double T;
    double B;
    unsigned long long step;
    double energy;
    double magnetisation;
} IsingLattice;

// Prototype functions declarations
IsingLattice *create_lattice(int N, double T, double B);
void initialize_lattice(IsingLattice *lattice);
double ising_hamiltonian(IsingLattice *lattice);
double ising_magnetisation(IsingLattice *lattice);
double delta_energy(IsingLattice *lattice, int i, int j);
int metropolis_algorithm(IsingLattice *lattice);
void write_csv(IsingLattice *lattice);

int main() {
    srand((unsigned)time(NULL)); // Random number generator
    int CASE; // Variable to hold user choice
    int N;
    double T;
    double B;
    char proceed;

    printf("2D Ising Model Simulation\n");
    printf("-------------------------\n");
    printf("Select N dimension for square lattice: \n");
    printf("1) 512x512  2) 256x256  3) Custom (Enter digit of choice)\n>");
    if (scanf("%d", &CASE)!=1){
        printf("Error: Invalid Input \n");
        return 1;
    }
    switch (CASE){ // Switch case for user input of lattice size
        case 1:
            N = 512;
            break;
        case 2:
            N=256;
            break;
        case 3:
            printf("Enter custom square lattice dimenion:\n>");
            if (scanf("%d", &N) != 1) return 1;
            break;
        default:
            printf("Error: invalid choice.\n");
            return 1;
    }
    printf("\nSelect temperature T for simulation:\n");
    printf("1) Low T  2) Critical T  3) High T  4) Custom\n>");
    if (scanf("%d", &CASE) != 1) return 1;
    switch (CASE){ // Switch case for user input of temperature
        case 1:
            T = 1.5;
            break;
        case 2:
            T = 2.269; // Critical temperature for 2D Ising model
            break;
        case 3:
            T = 7.5;
            break;
        case 4:
            printf("Enter custom temperature: (J/k unitless)\n>");
            if (scanf("%lf", &T) != 1) return 1;
            break;
        default:
            printf("Error: invalid choice.\n");
            return 1;
    }
    printf("\nSelect external magnetic field B:\n");
    printf("1) No field  2) Field on  3) Custom \n>");
    if (scanf("%d", &CASE) != 1) return 1;
    switch (CASE){ // Switch case for user input of external magnetic field
        case 1:
            B = 0;
            break;
        case 2:
            B = 0.1;
            break;
        case 3:
            printf("Enter custom field strength:\n>");
            if (scanf("%lf", &B) != 1) return 1;
            break;
        default:
            printf("Error: invalid choice.\n");
            return 1;
    }

    printf("Begin Ising simulation with N = %d, T = %.2f, B = %.2f? \n Y/N \n>", N, T, B); // Confirm user inputs before proceeding
    scanf(" %c", &proceed);

    if (proceed == 'Y' || proceed == 'y') { // Proceed with simulation if user confirms, allowing for case insensitivity
        IsingLattice *lattice = create_lattice(N, T, B); // Create lattice with user parameters
        initialize_lattice(lattice); // Initialize lattice with random spins

        int total_steps = 100000; // Total Monte Carlo steps for simulation
        for (int step = 0; step < total_steps; step++) {
            for (int attempt = 0; attempt < N * N; attempt++) { //Sweep through entire lattice per each monte carlo step
                metropolis_algorithm(lattice);
            }

            lattice->step = step; // Update current step count

            if (step % 200 == 0) { // Write to CSV every 200 steps meaning 500 frames for 100000 steps
                write_csv(lattice);
            }
        }

        printf("Simulation completed. Data written to ising_data.csv\n");

        free(lattice->spin); // Free allocated memory for spins
        free(lattice); // Free allocated memory for lattice structure
    } else {
        printf("Simulation aborted.\n");
    }

    return 0;
}

IsingLattice *create_lattice(int N, double T, double B) {
    /*
    * Function to create and allocate memory for an Ising lattice structure.
    * Initializes lattice parameters and returns a pointer to the created IsingLattice.
    */

    IsingLattice *lattice = malloc(sizeof(IsingLattice));
    if (lattice == NULL) { // Error handling for memory allocation, check if malloc was successful
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    lattice->N = N;
    lattice->T = T;
    lattice->B = B;
    lattice->step = 0;
    lattice->energy = 0.0;
    lattice->magnetisation = 0.0;
    lattice->spin = malloc(N * N * sizeof(Spin)); // Allocate memory for N x N spins with spin of char size
 
    if (lattice->spin == NULL) {
        fprintf(stderr, "Memory allocation for spin failed\n");
        free(lattice); // Free previously allocated memory before exiting
        exit(EXIT_FAILURE);
    }

    return lattice;
}

void initialize_lattice(IsingLattice *lattice) {
    /*
    * Function to initialize the lattice with random spin values (+1 or -1).
    * Also calculates initial energy and magnetisation.
    */

    int N = lattice->N;

    for (int i = 0; i < N * N; i++) {
        lattice->spin[i] = (rand() % 2) ? SPIN_UP : SPIN_DOWN; // Randomly assign spin up or down
    }

    lattice->energy = ising_hamiltonian(lattice); // Calculate initial energy
    lattice->magnetisation = ising_magnetisation(lattice); // Calculate initial magnetisation
}

double ising_hamiltonian(IsingLattice *lattice) {
    /*
    * Function to calculate the total energy of the lattice using the Ising Hamiltonian.
    * Considers nearest-neighbour interactions and external magnetic field.
    */

    int N = lattice->N;
    double energy = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { //Iterate through row and columns of spin lattice to calculate energy
            Spin s = lattice->spin[i * N + j];

            Spin right = lattice->spin[i * N + ((j + 1) % N)];
            Spin down  = lattice->spin[((i + 1) % N) * N + j];
            energy -= s * (right + down); // Nearest-neighbour interaction contribution
            energy -= lattice->B * s; // External magnetic field contribution
        }
    }
    return energy;
}

double ising_magnetisation(IsingLattice *lattice) {
    /*
    * Function to calculate the normalised magnetisation of the lattice.
    * Returns the average magnetisation per spin.
    */

    int N = lattice->N;
    double total_magnetisation = 0.0;

    for (int i = 0; i < N * N; i++) { 
        total_magnetisation += lattice->spin[i]; // Sum all spin values
    }

    return total_magnetisation / (N * N); // Normalise by total number of spins
}

double delta_energy(IsingLattice *lattice, int i, int j) {
    /*
    * Function to calculate the change in energy (ΔE) if the spin at (i, j) is flipped.
    * Takes i and j as arguments to identify the spin position in the lattice.
    */

    int N = lattice->N;
    Spin spin = lattice->spin[i * N + j];

    Spin right = lattice->spin[i * N + ((j + 1) % N)];
    Spin left  = lattice->spin[i * N + ((j - 1 + N) % N)];
    Spin down  = lattice->spin[((i + 1) % N) * N + j];
    Spin up    = lattice->spin[((i - 1 + N) % N) * N + j]; 

    double delta = 2.0 * spin * (right + left + up + down + lattice->B); // Calculate energy change for flipping the spin
    return delta;
}

int metropolis_algorithm(IsingLattice *lattice) {
    /*
    * Function to perform a single Metropolis algorithm step.
    * Randomly selects a spin, calculates the energy change if flipped,
    * and decides whether to flip the spin based on the Metropolis criteria.
    */

    int N = lattice->N;
    int i = rand() % N; // Randomly select spin positions whose value is less than N
    int j = rand() % N;

    double dE = delta_energy(lattice, i, j);

    if (dE <= 0.0) { // Energy change is always favourable when dE <= 0
        lattice->spin[i * N + j] *= -1; // Flip spin state
        return 1;
    } else {
        if (lattice->T <= 0.0) return 0;

        double acceptance_prob = exp(-dE / lattice->T); //Metroplis criteria for unfavourable energy changes
        double r = (double)rand() / (double)RAND_MAX; // Random probability value between 0 and 1
        if (r < acceptance_prob) {
            lattice->spin[i * N + j] *= -1;
            return 1;
        }
    }
    return 0;
}

void write_csv(IsingLattice *lattice) {
    /*
    * Function to write the current state of the lattice to a CSV file.
    * Each row contains: step, energy, magnetisation, followed by the spin states (0 or 1).
    */

    FILE *fp = fopen("ising_data.csv", "a");
    if (fp == NULL) { // Error handling for file opening, check if file opened successfully
        perror("Error opening CSV"); 
        exit(EXIT_FAILURE);
    }

    int N = lattice->N;

    lattice->energy = ising_hamiltonian(lattice);
    lattice->magnetisation = ising_magnetisation(lattice);

    fprintf(fp, "%llu,%f,%f",
            lattice->step,
            lattice->energy,
            lattice->magnetisation);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { //Iterate through row and columns to write spins states

            Spin s = lattice->spin[i * N + j];
            int bit = (s + 1) / 2; // Convert spin to 0 or 1 for CSV binary representation

            fprintf(fp, ",%d", bit);
        }
    }

    fprintf(fp, "\n");
    fclose(fp);
}

