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
double delta_energy(IsingLattice *lattice, int i, int j);
double metropolis_algorithm(IsingLattice *lattice);


int main(){
    srand(time(NULL));
    int N;
    double T;
    double B;
    char proceed;

    printf("2D Ising Model Simulation\n");
    printf("-------------------------\n");
    printf("Enter N dimension for square lattice: \n");
    scanf("%d", &N);
    printf("Enter Temperature T: \n");
    scanf("%lf", &T);
    printf("Enter External Magnetic Field B: \n");
    scanf("%lf", &B);
    printf("Begin simulation with N=%d, T=%f, B=%f? \n Y/N \n", N, T, B);
    scanf(" %c", &proceed);

    switch (proceed){
        case 'Y':{
             IsingLattice *lattice = create_lattice(N, T, B);

            double initial_energy = ising_hamiltonian(lattice);
            double initial_magnetisation = ising_magnetisation(lattice);
            printf("Initial Energy: %f\n", initial_energy);
            printf("Initial Magnetisation: %f\n", initial_magnetisation);

            for (int step = 0; step < 100000; step++){
                metropolis_algorithm(lattice);
            }

            double final_energy = ising_hamiltonian(lattice);
            double final_magnetisation = ising_magnetisation(lattice);
            printf("Final Energy: %f\n", final_energy);
            printf("Final Magnetisation: %f\n", final_magnetisation);


            break;
        }
        case 'N':{
            printf("Simulation aborted.\n");
            break;
        }
        default:{
            printf("ERROR: Invalid option\n");
            break;
        }
    }
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

double delta_energy(IsingLattice *lattice, int i, int j){
    int N = lattice->N;
    double B = lattice->B;

    Spin spin = lattice->spin[i * N + j];
    Spin right_adj = lattice->spin[i * N + ((j + 1) % N)];
    Spin left_adj  = lattice->spin[i * N + ((j - 1 + N) % N)];
    Spin down_adj  = lattice->spin[((i + 1) % N) * N + j];
    Spin up_adj    = lattice->spin[((i - 1 + N) % N) * N + j];

    double delta_energy = 2 * spin * (right_adj + left_adj + down_adj + up_adj +  B);

    return delta_energy;
}

double metropolis_algorithm(IsingLattice *lattice){
    int N = lattice->N;
    double T = lattice->T;

    int i = rand() % N;
    int j = rand() % N;

    double energy = delta_energy(lattice, i, j);

    if (energy <= 0){
        lattice->spin[i * N + j] *= -1;
        lattice->step += 1;
        return energy;
    } else {
        double acceptance_prob = exp(-energy / T);
        double r = (double)rand() / RAND_MAX;
        if (r < acceptance_prob){
            lattice->spin[i * N + j] *= -1;
            lattice->step += 1;
            return energy;
        } else {
            return 0.0; 
        }
    }
}