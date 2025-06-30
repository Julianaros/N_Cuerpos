#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <mpi.h>

typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
    float ax, ay, az;  // Current accelerations
    float ax_prev, ay_prev, az_prev;  // Previous accelerations for Verlet
} Particle;

/* Function Declarations */
int convertStringToInt(char *str);
void randomizeParticles(Particle *particles, int n, int rank, int size);

int main(int argc, char *argv[]) {
    int num_particles = 100000; // Default number of particles if no parameter is provided on the command line
    int rank, size;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc > 1) {
        // If a command-line parameter indicating the number of particles is provided, convert it to an integer
        num_particles = convertStringToInt(argv[1]);
    }
    
    if (rank == 0) {
        printf("MPI particle generation using %d processes\n", size);
    }
    
    // Allocate memory for particles
    Particle *particles = NULL;
    particles = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        printf("Process %d: Error: Unable to allocate memory for particles\n", rank);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    // Start timer
    double start_time = MPI_Wtime();
    
    // Randomize particles with MPI parallelization
    randomizeParticles(particles, num_particles, rank, size);
    
    // Stop timer
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    
    // Only rank 0 writes to file
    if (rank == 0) {
        // Write particles to file in binary mode
        FILE *file = fopen("particles.txt", "wb");
        if (file != NULL) {
            size_t written = fwrite(particles, sizeof(Particle), num_particles, file);
            fclose(file);
            
            if (written != num_particles) {
                printf("Error: Could not write all particles to file\n");
                free(particles);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
            
            printf("%d particles have been created with random values and written to file: particles.txt in binary format.\n", num_particles);
            printf("Elapsed time for generating %d particles is: %f seconds.\n", num_particles, elapsed_time);
            printf("Particle structure includes accelerations for Verlet integration.\n");
            printf("Processes used: %d\n", size);
        } else {
            printf("Error: Unable to open particles.txt for writing\n");
            free(particles);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
    
    // Free allocated memory
    free(particles);
    MPI_Finalize();
    return 0;
}

/* Convert string to integer */
int convertStringToInt(char *str) {
    char *endptr;
    long val;
    errno = 0; // To distinguish success/failure after the call
    
    val = strtol(str, &endptr, 10);
    
    /* Check for possible errors */
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0)) {
        perror("strtol");
        exit(EXIT_FAILURE);
    }
    
    if (endptr == str) {
        fprintf(stderr, "No digits were found\n");
        exit(EXIT_FAILURE);
    }
    
    /* If we are here, strtol() successfully converted a number */
    return (int)val;
}

/* Initialize particle state with random values using MPI */
void randomizeParticles(Particle *particles, int n, int rank, int size) {
    int particles_per_process = n / size;
    int remaining_particles = n % size;
    
    // Calculate start and end indices for this process
    int start = rank * particles_per_process;
    int end = (rank + 1) * particles_per_process;
    
    // Distribute remaining particles to the last processes
    if (rank < remaining_particles) {
        start += rank;
        end += rank + 1;
    } else {
        start += remaining_particles;
        end += remaining_particles;
    }
    
    // Seed random number generator differently for each process
    srand(time(NULL) + rank);
    
    // Initialize particles assigned to this process
    for (int i = start; i < end; i++) {
        particles[i].mass = 2.0f; // Arbitrarily chosen value for particle mass -> 2.0
        particles[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f; // Random number between -1 and 1
        particles[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        particles[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        particles[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        particles[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        particles[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        
        // Initialize accelerations to zero
        particles[i].ax = 0.0f;
        particles[i].ay = 0.0f;
        particles[i].az = 0.0f;
        particles[i].ax_prev = 0.0f;
        particles[i].ay_prev = 0.0f;
        particles[i].az_prev = 0.0f;
    }
    
    // Gather all particles to rank 0
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particles, 
                  particles_per_process * sizeof(Particle), MPI_BYTE, MPI_COMM_WORLD);
    
    // Handle remaining particles separately if needed
    if (remaining_particles > 0) {
        // Create a temporary array for remaining particles
        Particle *temp_particles = malloc(remaining_particles * sizeof(Particle));
        
        if (rank < remaining_particles) {
            // Copy the extra particle this process generated
            temp_particles[rank] = particles[end - 1];
        }
        
        // Gather remaining particles
        MPI_Allgather(temp_particles, sizeof(Particle), MPI_BYTE,
                      &particles[size * particles_per_process], sizeof(Particle), MPI_BYTE,
                      MPI_COMM_WORLD);
        
        free(temp_particles);
    }
}