#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>

#define SOFTENING 1e-9f

/* Pthreads parallel implementation of the simulation of the n-body problem
   Using Velocity Verlet integration method with pthreads parallelization
   Based on the sequential Verlet implementation but with pthreads optimizations */

typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
    float ax, ay, az;  // Current acceleration
    float ax_prev, ay_prev, az_prev;  // Previous acceleration for Verlet
} Particle;

typedef struct {
    Particle *particles;
    int start;
    int end;
    int n;
    float dt;
} ThreadData;

/* Global variables for thread synchronization */
static int num_threads;
static pthread_barrier_t barrier;

/* Function declarations */
int convertStringToInt(char *str);
void computeAccelerations(Particle *p, int n);
void verletIntegration(Particle *p, float dt, int n);
void* computeAccelerationsThread(void* arg);
void* verletIntegrationThread(void* arg);

int main(int argc, char* argv[]) {

    int nBodies = 1000; // Default number of bodies if no parameters are given from the command line
    if (argc > 1) nBodies = convertStringToInt(argv[1]);

    if (argc > 2) {
        num_threads = convertStringToInt(argv[2]);
    } else {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads <= 0) num_threads = 4;
    }

    const float dt = 0.01f; // Time step
    const int nIters = 10;  // Simulation iterations

    struct timespec startIter, endIter, startTotal, endTotal;
    double totalTime = 0.0;

    // Initialize barrier for thread synchronization
    if (pthread_barrier_init(&barrier, NULL, num_threads) != 0) {
        printf("Error initializing barrier\n");
        exit(EXIT_FAILURE);
    }

    // First generate particles using pthreads particle generator
    system("gcc -pthread -o particle_production_pthreads particle_production_pthreads.c -lm -O2");
    char command[256];
    sprintf(command, "./particle_production_pthreads %d %d", nBodies, num_threads);
    system(command);

    Particle *particles = NULL;
    particles = (Particle *)malloc(nBodies * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Unable to allocate memory for particles\n");
        exit(EXIT_FAILURE);
    }

    FILE *fileRead = fopen("particles.txt", "rb");
    if (fileRead == NULL) {
        /* Unable to open the file */
        printf("\nUnable to open the file particles.txt.\n");
        free(particles);
        exit(EXIT_FAILURE);
    }

    size_t particlesRead = fread(particles, sizeof(Particle), nBodies, fileRead);
    if (particlesRead != nBodies) {
        /* The number of particles to read is different than expected */
        printf("ERROR: Expected to read %d particles, but read %zu particles from file\n", nBodies, particlesRead);
        fclose(fileRead);
        free(particles);
        exit(EXIT_FAILURE);
    }
    fclose(fileRead);

    // Initialize accelerations to zero
    for (int i = 0; i < nBodies; i++) {
        particles[i].ax = particles[i].ay = particles[i].az = 0.0f;
        particles[i].ax_prev = particles[i].ay_prev = particles[i].az_prev = 0.0f;
    }
    
    // Compute initial accelerations
    computeAccelerations(particles, nBodies);

    printf("\n=== N-Body Simulation Results (Pthreads Velocity Verlet) ===\n");
    printf("Number of particles: %d\n", nBodies);
    printf("Time step: %f\n", dt);
    printf("Number of iterations: %d\n", nIters);
    printf("Integration method: Velocity Verlet (Pthreads)\n");
    printf("Threads used: %d\n", num_threads);
    printf("===========================================================\n");

    clock_gettime(CLOCK_MONOTONIC, &startTotal);

    for (int iter = 1; iter <= nIters; iter++) {
        clock_gettime(CLOCK_MONOTONIC, &startIter);

        verletIntegration(particles, dt, nBodies); // Pthreads Velocity Verlet integration

        clock_gettime(CLOCK_MONOTONIC, &endIter);
        double iterTime = (endIter.tv_sec - startIter.tv_sec) + (endIter.tv_nsec - startIter.tv_nsec) / 1e9;
        printf("Iteration %d of %d completed in %f seconds\n", iter, nIters, iterTime);
    }

    clock_gettime(CLOCK_MONOTONIC, &endTotal);
    totalTime = (endTotal.tv_sec - startTotal.tv_sec) + (endTotal.tv_nsec - startTotal.tv_nsec) / 1e9;
    double avgTime = totalTime / (double)(nIters);
    
    printf("\n=== Performance Summary ===\n");
    printf("Avg iteration time: %f seconds\n", avgTime);
    printf("Total simulation time: %f seconds\n", totalTime);
    printf("Particles processed: %d\n", nBodies);
    printf("Integration method: Velocity Verlet (Pthreads)\n");
    printf("Threads used: %d\n", num_threads);
    printf("==========================\n\n");

    /* Write the output to a file to evaluate correctness by comparing with sequential output */
    FILE *fileWrite = fopen("parallel_output.txt", "wb");
    if (fileWrite != NULL) {
        size_t written = fwrite(particles, sizeof(Particle), nBodies, fileWrite);
        fclose(fileWrite);
        if (written == nBodies) {
            printf("Final particle states written to parallel_output.txt\n");
        } else {
            printf("Warning: Could not write all particles to output file\n");
        }
    } else {
        printf("Warning: Could not create output file parallel_output.txt\n");
    }

    // Clean up
    pthread_barrier_destroy(&barrier);
    free(particles);
    return 0;
}

/* Conversion from string to integer */
int convertStringToInt(char *str) {
    char *endptr;
    long val;  
    errno = 0;  // To distinguish success/failure after the call

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

    /* If we are here, strtol() has converted a number correctly */
    return (int)val;
}

/* Thread function for computing accelerations */
void* computeAccelerationsThread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    Particle* p = data->particles;
    int n = data->n;
    
    // Reset accelerations for assigned particles
    for (int i = data->start; i < data->end; i++) {
        p[i].ax = p[i].ay = p[i].az = 0.0f;
    }
    
    // Synchronize all threads before computing forces
    pthread_barrier_wait(&barrier);
    
    // Compute gravitational forces and accelerations
    for (int i = data->start; i < data->end; i++) {
        float ax_local = 0.0f, ay_local = 0.0f, az_local = 0.0f;
        
        for (int j = 0; j < n; j++) {
            if (i != j) { // Don't compute force on self
                float dx = p[j].x - p[i].x;
                float dy = p[j].y - p[i].y;
                float dz = p[j].z - p[i].z;
                float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
                float invDist = 1.0f / sqrtf(distSqr);
                float invDist3 = invDist * invDist * invDist;

                // F = G * m1 * m2 / r^2, but we normalize G=1 for simplicity
                // a = F/m = G * m_other / r^2
                float acceleration = p[j].mass * invDist3;
                
                ax_local += acceleration * dx;
                ay_local += acceleration * dy;
                az_local += acceleration * dz;
            }
        }
        
        // Store the computed accelerations
        p[i].ax = ax_local;
        p[i].ay = ay_local;
        p[i].az = az_local;
    }
    
    return NULL;
}

/* Function that computes accelerations for all particles using pthreads */
void computeAccelerations(Particle *p, int n) {
    pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
    ThreadData* thread_data = malloc(num_threads * sizeof(ThreadData));
    
    int particles_per_thread = n / num_threads;
    int remaining_particles = n % num_threads;
    
    // Create threads
    for (int t = 0; t < num_threads; t++) {
        thread_data[t].particles = p;
        thread_data[t].n = n;
        thread_data[t].start = t * particles_per_thread;
        thread_data[t].end = (t + 1) * particles_per_thread;
        
        // Distribute remaining particles to last thread
        if (t == num_threads - 1) {
            thread_data[t].end += remaining_particles;
        }
        
        if (pthread_create(&threads[t], NULL, computeAccelerationsThread, &thread_data[t]) != 0) {
            printf("Error creating thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    // Wait for all threads to complete
    for (int t = 0; t < num_threads; t++) {
        if (pthread_join(threads[t], NULL) != 0) {
            printf("Error joining thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    free(threads);
    free(thread_data);
}

/* Thread function for Verlet integration */
void* verletIntegrationThread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    Particle* p = data->particles;
    float dt = data->dt;
    int n = data->n;
    
    // Store current accelerations as previous
    for (int i = data->start; i < data->end; i++) {
        p[i].ax_prev = p[i].ax;
        p[i].ay_prev = p[i].ay;
        p[i].az_prev = p[i].az;
    }
    
    // Synchronize all threads
    pthread_barrier_wait(&barrier);
    
    // Update positions: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    for (int i = data->start; i < data->end; i++) {
        p[i].x += p[i].vx * dt + 0.5f * p[i].ax * dt * dt;
        p[i].y += p[i].vy * dt + 0.5f * p[i].ay * dt * dt;
        p[i].z += p[i].vz * dt + 0.5f * p[i].az * dt * dt;
    }
    
    // Synchronize before computing new accelerations
    pthread_barrier_wait(&barrier);
    
    return NULL;
}

/* Velocity Verlet integration step with pthreads parallelization */
void verletIntegration(Particle *p, float dt, int n) {
    pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
    ThreadData* thread_data = malloc(num_threads * sizeof(ThreadData));
    
    int particles_per_thread = n / num_threads;
    int remaining_particles = n % num_threads;
    
    // Create threads for position update
    for (int t = 0; t < num_threads; t++) {
        thread_data[t].particles = p;
        thread_data[t].n = n;
        thread_data[t].dt = dt;
        thread_data[t].start = t * particles_per_thread;
        thread_data[t].end = (t + 1) * particles_per_thread;
        
        // Distribute remaining particles to last thread
        if (t == num_threads - 1) {
            thread_data[t].end += remaining_particles;
        }
        
        if (pthread_create(&threads[t], NULL, verletIntegrationThread, &thread_data[t]) != 0) {
            printf("Error creating thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    // Wait for position update to complete
    for (int t = 0; t < num_threads; t++) {
        if (pthread_join(threads[t], NULL) != 0) {
            printf("Error joining thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    // Compute new accelerations at t+dt
    computeAccelerations(p, n);
    
    // Update velocities: v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]*dt
    for (int t = 0; t < num_threads; t++) {
        if (pthread_create(&threads[t], NULL, verletIntegrationThread, &thread_data[t]) != 0) {
            printf("Error creating thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    // Custom thread function for velocity update
    for (int t = 0; t < num_threads; t++) {
        if (pthread_join(threads[t], NULL) != 0) {
            printf("Error joining thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }
    
    // Update velocities in main thread (simpler approach)
    for (int i = 0; i < n; i++) {
        p[i].vx += 0.5f * (p[i].ax_prev + p[i].ax) * dt;
        p[i].vy += 0.5f * (p[i].ay_prev + p[i].ay) * dt;
        p[i].vz += 0.5f * (p[i].az_prev + p[i].az) * dt;
    }
    
    free(threads);
    free(thread_data);
}