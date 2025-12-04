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
    int N;
    double T;
    double B;
    char proceed;

    printf("2D Ising Model Simulation\n");
    printf("-------------------------\n");
    printf("Enter N dimension for square lattice: \n");
    if (scanf("%d", &N) != 1) return 1;
    printf("Enter Temperature T: \n");
    if (scanf("%lf", &T) != 1) return 1;
    printf("Enter External Magnetic Field B: \n");
    if (scanf("%lf", &B) != 1) return 1;
    printf("Begin simulation with N=%d, T=%f, B=%f? \n Y/N \n", N, T, B);
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

