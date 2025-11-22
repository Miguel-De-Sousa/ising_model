#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef char Spin;
static const Spin SPIN_UP   =  1;
static const Spin SPIN_DOWN = -1;

typedef struct{
    int N;
    Spin *spin;
    double T;
    double B;
    unsigned long long step;
    double energy;
    double magnetisation;
} IsingLattice;

IsingLattice *create_lattice(int N, double T, double B);
double ising_hamiltonian(IsingLattice *lattice);
double ising_magnetisation(IsingLattice *lattice);

int main(){
    srand(time(NULL));

    IsingLattice *lattice = create_lattice(20, 2.0, 0.5);
    double E = ising_hamiltonian(lattice);
    double M = ising_magnetisation(lattice);

    printf("Total Energy = %fJ\n", E);
    printf("Magnetisation = %fT\n", M);

    return 0;
}

IsingLattice *create_lattice(int N, double T, double B){

    IsingLattice *lattice = malloc(sizeof(IsingLattice));

    if (lattice == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    lattice->N = N;
    lattice->T = T;
    lattice->B = B;
    lattice->step = 0;
    lattice->spin = malloc(N * N * sizeof(char));

    if (lattice->spin == NULL){
        fprintf(stderr, "Memory allocation for spin failed\n");
        free(lattice);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N * N; i++){
        lattice->spin[i] = (rand() % 2) ? SPIN_UP : SPIN_DOWN;
    }

    return lattice;
}

double ising_hamiltonian(IsingLattice *lattice){
    int N = lattice->N;
    double B = lattice->B;

    double adj_energy = 0;
    double field_energy = 0;

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            Spin spin = lattice->spin[i * N + j];

            Spin right_adj = lattice->spin[i * N + ((j + 1) % N)];
            Spin down_adj = lattice->spin[((i + 1) % N) * N + j];

            adj_energy -= spin * (right_adj + down_adj);
            field_energy -= B * spin;
        }
    }
    double total_energy = adj_energy + field_energy;
    lattice->energy = total_energy;
    return total_energy;
}

double ising_magnetisation(IsingLattice *lattice){
    int N = lattice->N;
    double total_magnetisation = 0;

    for (int i=0; i<N*N; i++){
        total_magnetisation += lattice->spin[i];
    }

    double normalised_mag = total_magnetisation / (N * N);
    lattice->magnetisation = normalised_mag;
    return normalised_mag;
}
