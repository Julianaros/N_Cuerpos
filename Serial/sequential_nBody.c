#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#define SOFTENING 1e-9f

/* Implementation of the simulation of the n-body problem
   Sequential version, taken from the example given 
   at the link https://github.com/harrism/mini-nbody/blob/master/nbody.c
   and adapted to the Linux environment */

typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;

/* Function declarations */
int convertStringToInt(char *str);
void bodyForce(Particle *p, float dt, int n);

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

    printf("\n=== N-Body Simulation Results ===\n");
    printf("Number of particles: %d\n", nBodies);
    printf("Time step: %f\n", dt);
    printf("Number of iterations: %d\n", nIters);
    printf("==================================\n");

    for (int iter = 1; iter <= nIters; iter++) {
        startIter = clock();

        bodyForce(particles, dt, nBodies); // Compute inter-body forces

        for (int i = 0; i < nBodies; i++) { // Integrate position
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;
            particles[i].z += particles[i].vz * dt;
        }

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
    printf("==========================\n\n");

    /* Write the output to a file to evaluate correctness by comparing with parallel output */
    FILE *fileWrite = fopen("sequential_output.txt", "wb");
    if (fileWrite != NULL) {
        size_t written = fwrite(particles, sizeof(Particle), nBodies, fileWrite);
        fclose(fileWrite);
        if (written == nBodies) {
            printf("Final particle states written to sequential_output.txt\n");
        } else {
            printf("Warning: Could not write all particles to output file\n");
        }
    } else {
        printf("Warning: Could not create output file sequential_output.txt\n");
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

/* Function that performs computation */
void bodyForce(Particle *p, float dt, int n) {
    for (int i = 0; i < n; i++) {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < n; j++) {
            if (i != j) { // Don't compute force on self
                float dx = p[j].x - p[i].x;
                float dy = p[j].y - p[i].y;
                float dz = p[j].z - p[i].z;
                float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
                float invDist = 1.0f / sqrtf(distSqr);
                float invDist3 = invDist * invDist * invDist;

                Fx += p[j].mass * dx * invDist3;
                Fy += p[j].mass * dy * invDist3;
                Fz += p[j].mass * dz * invDist3;
            }
        }

        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }
}