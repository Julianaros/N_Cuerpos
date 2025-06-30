#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <omp.h>

typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
    float ax, ay, az;  // Current accelerations
    float ax_prev, ay_prev, az_prev;  // Previous accelerations for Verlet
} Particle;

/* Function Declarations */
int convertStringToInt(char *str);
void randomizeParticles(Particle *particles, int n);

int main(int argc, char *argv[]) {
    int num_particles = 100000; // Default number of particles if no parameter is provided on the command line
    
    if (argc > 1) {
        // If a command-line parameter indicating the number of particles is provided, convert it to an integer
        num_particles = convertStringToInt(argv[1]);
    }
    
    // Get number of available threads
    int num_threads = omp_get_max_threads();
    printf("OpenMP particle generation using %d threads\n", num_threads);
    
    // Allocate memory for particles
    Particle *particles = NULL;
    particles = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Unable to allocate memory for particles\n");
        exit(EXIT_FAILURE);
    }
    
    // Start timer
    double start = omp_get_wtime();
    
    // Randomize particles with OpenMP parallelization
    randomizeParticles(particles, num_particles);
    
    // Write particles to file in binary mode
    FILE *file = fopen("particles.txt", "wb");
    if (file != NULL) {
        size_t written = fwrite(particles, sizeof(Particle), num_particles, file);
        fclose(file);
        
        if (written != num_particles) {
            printf("Error: Could not write all particles to file\n");
            free(particles);
            exit(EXIT_FAILURE);
        }
        
        // Stop timer
        double end = omp_get_wtime();
        double elapsed_time = end - start;
        
        printf("%d particles have been created with random values and written to file: particles.txt in binary format.\n", num_particles);
        printf("Elapsed time for generating %d particles is: %f seconds.\n", num_particles, elapsed_time);
        printf("Particle structure includes accelerations for Verlet integration.\n");
        printf("Threads used: %d\n", num_threads);
    } else {
        printf("Error: Unable to open particles.txt for writing\n");
        free(particles);
        exit(EXIT_FAILURE);
    }
    
    // Free allocated memory
    free(particles);
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

/* Initialize particle state with random values using OpenMP */
void randomizeParticles(Particle *particles, int n) {
    // Use OpenMP to parallelize particle initialization
    #pragma omp parallel
    {
        // Each thread gets its own random seed based on thread ID
        unsigned int seed = omp_get_thread_num();
        
        #pragma omp for
        for (int i = 0; i < n; i++) {
            particles[i].mass = 2.0f; // Arbitrarily chosen value for particle mass -> 2.0
            particles[i].x = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f; // Random number between -1 and 1
            particles[i].y = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
            particles[i].z = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
            particles[i].vx = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
            particles[i].vy = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
            particles[i].vz = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
            
            // Initialize accelerations to zero
            particles[i].ax = 0.0f;
            particles[i].ay = 0.0f;
            particles[i].az = 0.0f;
            particles[i].ax_prev = 0.0f;
            particles[i].ay_prev = 0.0f;
            particles[i].az_prev = 0.0f;
        }
    }
}