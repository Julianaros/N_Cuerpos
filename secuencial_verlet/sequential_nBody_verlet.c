#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#define SOFTENING 1e-9f

/* Implementation of the simulation of the n-body problem
   Sequential version using Velocity Verlet integration method
   Based on the original Euler implementation but with improved
   numerical stability and accuracy */

typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
    float ax, ay, az;  // Current acceleration
    float ax_prev, ay_prev, az_prev;  // Previous acceleration for Verlet
} Particle;

/* Function declarations */
int convertStringToInt(char *str);
void computeAccelerations(Particle *p, int n);
void verletIntegration(Particle *p, float dt, int n);

int main(int argc, char* argv[]) {

    int nBodies = 1000; // Default number of bodies if no parameters are given from the command line
    if (argc > 1) nBodies = convertStringToInt(argv[1]);

    const float dt = 0.01f; // Time step
    const int nIters = 10;  // Simulation iterations

    clock_t startIter, endIter;
    clock_t startTotal = clock(), endTotal;
    double totalTime = 0.0;

    // First generate particles
    system("gcc -o particle_production particle_production.c -lm");
    char command[256];
    sprintf(command, "./particle_production %d", nBodies);
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

    // Initialize accelerations to zero and compute initial accelerations
    for (int i = 0; i < nBodies; i++) {
        particles[i].ax = particles[i].ay = particles[i].az = 0.0f;
        particles[i].ax_prev = particles[i].ay_prev = particles[i].az_prev = 0.0f;
    }
    
    // Compute initial accelerations
    computeAccelerations(particles, nBodies);

    printf("\n=== N-Body Simulation Results (Velocity Verlet) ===\n");
    printf("Number of particles: %d\n", nBodies);
    printf("Time step: %f\n", dt);
    printf("Number of iterations: %d\n", nIters);
    printf("Integration method: Velocity Verlet\n");
    printf("===================================================\n");

    for (int iter = 1; iter <= nIters; iter++) {
        startIter = clock();

        verletIntegration(particles, dt, nBodies); // Velocity Verlet integration

        endIter = clock() - startIter;
        printf("Iteration %d of %d completed in %f seconds\n", iter, nIters, (double)endIter / CLOCKS_PER_SEC);
    }

    endTotal = clock();
    totalTime = (double)(endTotal - startTotal) / CLOCKS_PER_SEC;
    double avgTime = totalTime / (double)(nIters);
    
    printf("\n=== Performance Summary ===\n");
    printf("Avg iteration time: %f seconds\n", avgTime);
    printf("Total simulation time: %f seconds\n", totalTime);
    printf("Particles processed: %d\n", nBodies);
    printf("Integration method: Velocity Verlet\n");
    printf("==========================\n\n");

    /* Write the output to a file to evaluate correctness by comparing with parallel output */
    FILE *fileWrite = fopen("sequential_verlet_output.txt", "wb");
    if (fileWrite != NULL) {
        size_t written = fwrite(particles, sizeof(Particle), nBodies, fileWrite);
        fclose(fileWrite);
        if (written == nBodies) {
            printf("Final particle states written to sequential_verlet_output.txt\n");
        } else {
            printf("Warning: Could not write all particles to output file\n");
        }
    } else {
        printf("Warning: Could not create output file sequential_verlet_output.txt\n");
    }

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

/* Function that computes accelerations for all particles */
void computeAccelerations(Particle *p, int n) {
    // Reset accelerations
    for (int i = 0; i < n; i++) {
        p[i].ax = p[i].ay = p[i].az = 0.0f;
    }
    
    // Compute gravitational forces and accelerations
    for (int i = 0; i < n; i++) {
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
                
                p[i].ax += acceleration * dx;
                p[i].ay += acceleration * dy;
                p[i].az += acceleration * dz;
            }
        }
    }
}

/* Velocity Verlet integration step */
void verletIntegration(Particle *p, float dt, int n) {
    // Store current accelerations as previous
    for (int i = 0; i < n; i++) {
        p[i].ax_prev = p[i].ax;
        p[i].ay_prev = p[i].ay;
        p[i].az_prev = p[i].az;
    }
    
    // Update positions: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    for (int i = 0; i < n; i++) {
        p[i].x += p[i].vx * dt + 0.5f * p[i].ax * dt * dt;
        p[i].y += p[i].vy * dt + 0.5f * p[i].ay * dt * dt;
        p[i].z += p[i].vz * dt + 0.5f * p[i].az * dt * dt;
    }
    
    // Compute new accelerations at t+dt
    computeAccelerations(p, n);
    
    // Update velocities: v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]*dt
    for (int i = 0; i < n; i++) {
        p[i].vx += 0.5f * (p[i].ax_prev + p[i].ax) * dt;
        p[i].vy += 0.5f * (p[i].ay_prev + p[i].ay) * dt;
        p[i].vz += 0.5f * (p[i].az_prev + p[i].az) * dt;
    }
}