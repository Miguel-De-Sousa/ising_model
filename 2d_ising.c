#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SPIN_UP 1
#define SPIN_DOWN -1

typedef struct{
    int N;
    short *spins;
    double T;
    unsigned long long step;
} square_Lattice;

square_Lattice * create_lattice(int N, double T){

    square_Lattice * lattice_malloc = malloc(sizeof(square_Lattice));

    if (lattice_malloc == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    lattice_malloc->N = N;
    lattice_malloc->T = T;
    lattice_malloc->step = 0;
    lattice_malloc->spins = malloc(N * N * sizeof(short));

    if (lattice_malloc->spins == NULL){
        fprintf(stderr, "Memory allocation for spins failed\n");
        free(lattice_malloc);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N * N; i++){
        lattice_malloc->spins[i] = (rand() % 2) ? SPIN_UP : SPIN_DOWN;
    }

    return lattice_malloc;
}