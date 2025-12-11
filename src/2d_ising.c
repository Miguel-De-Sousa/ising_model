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

typedef char Spin;
static const Spin SPIN_UP   =  1;
static const Spin SPIN_DOWN = -1;

typedef struct {
    int N;
    Spin *spin;
    double T;
    double B;
    unsigned long long step;
    double energy;
    double magnetisation;
} IsingLattice;

IsingLattice *create_lattice(int N, double T, double B);
void initialize_lattice(IsingLattice *lattice);
double ising_hamiltonian(IsingLattice *lattice);
double ising_magnetisation(IsingLattice *lattice);
double delta_energy(IsingLattice *lattice, int i, int j);
int metropolis_algorithm(IsingLattice *lattice);
void write_csv(IsingLattice *lattice);

int main() {
    srand((unsigned)time(NULL));
    int CASE;
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
    switch (CASE){
        case 1:
            N = 512;
            break;
        case 2:
            N=256;
            break;
        case 3:
            printf("Enter custom square lattice dimenion:\n>");
            if (scanf("%d", &CASE) != 1) return 1;
            N = CASE;
            break;
        default:
            printf("Error: invalid choice.\n");
            return 1;
    }
    printf("\nSelect temperature T for simulation:\n");
    printf("1) Low T  2) Critical T  3) High T  4) Custom\n>");
    if (scanf("%d", &CASE) != 1) return 1;
    switch (CASE){
        case 1:
            T = 1.5;
            break;
        case 2:
            T = 2.269;
            break;
        case 3:
            T = 7.5;
            break;
        case 4:
            printf("Enter custom temperature: (J/k unitless)\n>");
            if (scanf("%d", &CASE) != 1) return 1;
            T = CASE;
            break;
        default:
            printf("Error: invalid choice.\n");
            return 1;
    }
    printf("\nSelect external magnetic field B:\n");
    printf("1) No field  2) Field on  3) Custom \n>");
    if (scanf("%d", &CASE) != 1) return 1;
    switch (CASE){
        case 1:
            B = 0;
            break;
        case 2:
            B = 0.1;
            break;
        case 3:
            printf("Enter custom field strength:\n>");
            if (scanf("%d", &CASE) != 1) return 1;
            B = CASE;
            break;
        default:
            printf("Error: invalid choice.\n");
            return 1;
    }

    printf("Begin Ising simulation with N = %d, T = %.2f, B = %.2f? \n Y/N \n>", N, T, B);
    scanf(" %c", &proceed);

    if (proceed == 'Y' || proceed == 'y') {
        IsingLattice *lattice = create_lattice(N, T, B);
        initialize_lattice(lattice);

        int total_steps = 100000;
        for (int step = 0; step < total_steps; step++) {
            for (int attempt = 0; attempt < N * N; attempt++) {
                metropolis_algorithm(lattice);
            }

            lattice->step = step;

            if (step % 200 == 0) { 
                write_csv(lattice);
            }
        }

        printf("Simulation completed. Data written to ising_data.csv\n");

        free(lattice->spin);
        free(lattice);
    } else {
        printf("Simulation aborted.\n");
    }

    return 0;
}

IsingLattice *create_lattice(int N, double T, double B) {
    IsingLattice *lattice = malloc(sizeof(IsingLattice));
    if (lattice == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    lattice->N = N;
    lattice->T = T;
    lattice->B = B;
    lattice->step = 0;
    lattice->energy = 0.0;
    lattice->magnetisation = 0.0;
    lattice->spin = malloc(N * N * sizeof(Spin));

    if (lattice->spin == NULL) {
        fprintf(stderr, "Memory allocation for spin failed\n");
        free(lattice);
        exit(EXIT_FAILURE);
    }

    return lattice;
}

void initialize_lattice(IsingLattice *lattice) {
    int N = lattice->N;

    for (int i = 0; i < N * N; i++) {
        lattice->spin[i] = (rand() % 2) ? SPIN_UP : SPIN_DOWN;
    }

    lattice->energy = ising_hamiltonian(lattice);
    lattice->magnetisation = ising_magnetisation(lattice);
}

double ising_hamiltonian(IsingLattice *lattice) {
    int N = lattice->N;
    double energy = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Spin s = lattice->spin[i * N + j];

            Spin right = lattice->spin[i * N + ((j + 1) % N)];
            Spin down  = lattice->spin[((i + 1) % N) * N + j];
            energy -= s * (right + down);
            energy -= lattice->B * s;
        }
    }
    return energy;
}

double ising_magnetisation(IsingLattice *lattice) {
    int N = lattice->N;
    double total_magnetisation = 0.0;

    for (int i = 0; i < N * N; i++) {
        total_magnetisation += lattice->spin[i];
    }

    return total_magnetisation / (N * N);
}

double delta_energy(IsingLattice *lattice, int i, int j) {
    int N = lattice->N;
    Spin spin = lattice->spin[i * N + j];

    Spin right = lattice->spin[i * N + ((j + 1) % N)];
    Spin left  = lattice->spin[i * N + ((j - 1 + N) % N)];
    Spin down  = lattice->spin[((i + 1) % N) * N + j];
    Spin up    = lattice->spin[((i - 1 + N) % N) * N + j];

    double delta = 2.0 * spin * (right + left + up + down + lattice->B);
    return delta;
}

int metropolis_algorithm(IsingLattice *lattice) {
    int N = lattice->N;
    int i = rand() % N;
    int j = rand() % N;

    double dE = delta_energy(lattice, i, j);

    if (dE <= 0.0) {
        lattice->spin[i * N + j] *= -1;
        return 1;
    } else {
        if (lattice->T <= 0.0) return 0;

        double acceptance_prob = exp(-dE / lattice->T);
        double r = (double)rand() / (double)RAND_MAX;
        if (r < acceptance_prob) {
            lattice->spin[i * N + j] *= -1;
            return 1;
        }
    }
    return 0;
}

void write_csv(IsingLattice *lattice) {

    FILE *fp = fopen("ising_data.csv", "a");
    if (fp == NULL) {
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
        for (int j = 0; j < N; j++) {

            Spin s = lattice->spin[i * N + j];
            int bit = (s + 1) / 2;

            fprintf(fp, ",%d", bit);
        }
    }

    fprintf(fp, "\n");
    fclose(fp);
}

