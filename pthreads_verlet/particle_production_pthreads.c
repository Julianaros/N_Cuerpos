#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

// Global variable to store number of threads
int num_threads_global = 4;

typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
    float ax, ay, az;  // Current accelerations
    float ax_prev, ay_prev, az_prev;  // Previous accelerations for Verlet
} Particle;

typedef struct {
    Particle *particles;
    int start;
    int end;
    unsigned int seed;
} ThreadData;

/* Function Declarations */
int convertStringToInt(char *str);
void randomizeParticles(Particle *particles, int n);
void* randomizeParticlesThread(void* arg);

int main(int argc, char *argv[]) {
    int num_particles = 100000; // Default number of particles if no parameter is provided on the command line
    int num_threads = 4; // Default number of threads
    
    if (argc > 1) {
        // If a command-line parameter indicating the number of particles is provided, convert it to an integer
        num_particles = convertStringToInt(argv[1]);
    }
    
    if (argc > 2) {
        // If a second parameter is provided for number of threads
        num_threads = convertStringToInt(argv[2]);
    } else {
        // Get number of available threads (cores) if not specified
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads <= 0) num_threads = 4; // Default fallback
    }
    
    printf("Pthreads particle generation using %d threads\n", num_threads);
    
    // Set global thread count
    num_threads_global = num_threads;
    
    // Allocate memory for particles
    Particle *particles = NULL;
    particles = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Unable to allocate memory for particles\n");
        exit(EXIT_FAILURE);
    }
    
    // Start timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    // Randomize particles with pthreads parallelization
    randomizeParticles(particles, num_particles);
    
    // Stop timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
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

/* Thread function for initializing particles */
void* randomizeParticlesThread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    unsigned int seed = data->seed;
    
    for (int i = data->start; i < data->end; i++) {
        data->particles[i].mass = 2.0f; // Arbitrarily chosen value for particle mass -> 2.0
        data->particles[i].x = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f; // Random number between -1 and 1
        data->particles[i].y = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
        data->particles[i].z = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
        data->particles[i].vx = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
        data->particles[i].vy = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
        data->particles[i].vz = 2.0f * (rand_r(&seed) / (float)RAND_MAX) - 1.0f;
        
        // Initialize accelerations to zero
        data->particles[i].ax = 0.0f;
        data->particles[i].ay = 0.0f;
        data->particles[i].az = 0.0f;
        data->particles[i].ax_prev = 0.0f;
        data->particles[i].ay_prev = 0.0f;
        data->particles[i].az_prev = 0.0f;
    }
    
    return NULL;
}

/* Initialize particle state with random values using pthreads */
void randomizeParticles(Particle *particles, int n) {
    // Use global num_threads variable
    extern int num_threads_global;
    int threads_to_use = num_threads_global;
    
    pthread_t* threads = malloc(threads_to_use * sizeof(pthread_t));
    ThreadData* thread_data = malloc(threads_to_use * sizeof(ThreadData));
    
    int particles_per_thread = n / threads_to_use;
    int remaining_particles = n % threads_to_use;
    
    // Create threads
    for (int t = 0; t < threads_to_use; t++) {
        thread_data[t].particles = particles;
        thread_data[t].start = t * particles_per_thread;
        thread_data[t].end = (t + 1) * particles_per_thread;
        
        // Distribute remaining particles to last thread
        if (t == threads_to_use - 1) {
            thread_data[t].end += remaining_particles;
        }
        
        // Each thread gets a different seed
        thread_data[t].seed = time(NULL) + t;
        
        if (pthread_create(&threads[t], NULL, randomizeParticlesThread, &thread_data[t]) != 0) {
            printf("Error creating thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    // Wait for all threads to complete
    for (int t = 0; t < threads_to_use; t++) {
        if (pthread_join(threads[t], NULL) != 0) {
            printf("Error joining thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    free(threads);
    free(thread_data);
}