#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>   // For srand
#include <limits.h> // For LONG_MAX, LONG_MIN
#include <errno.h>  // For errno

#define SOFTENING 1e-9f

// Structure to represent a particle
typedef struct {
    float mass;
    float x, y, z;
    float vx, vy, vz;
    float ax, ay, az;        // Current accelerations
    float ax_prev, ay_prev, az_prev; // Previous accelerations for Verlet
} Particle;

// MPI custom datatype for Particle struct
MPI_Datatype MPI_Particle_Type;

/* Function Declarations */
int convertStringToInt(char *str);
void createMPITypeForParticle();
void randomizeLocalParticles(Particle *particles, int n_local, int rank);
void computeAccelerations(Particle *local_p, Particle *global_p, int local_n, int total_n, int start_idx);
void verletIntegration(Particle *local_p, Particle *global_p, float dt, int local_n, int total_n, int start_idx);

int main(int argc, char* argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Correctly get world_size

    int total_nBodies = 1000; // Default total number of bodies
    if (argc > 1) {
        total_nBodies = convertStringToInt(argv[1]);
    }

    const float dt = 0.01f; // Time step
    const int nIters = 10;  // Simulation iterations

    double start_total, end_total; // Total simulation time
    double start_iter, end_iter;   // Per-iteration time

    // Create custom MPI datatype for Particle struct
    createMPITypeForParticle();

    // Calculate local number of particles and their starting index
    int local_nBodies = total_nBodies / world_size;
    int start_idx = world_rank * local_nBodies;
    // Handle remainder particles for the last process
    if (world_rank == world_size - 1) {
        local_nBodies += (total_nBodies % world_size);
    }
    
    // Allocate memory for local particles
    Particle *local_particles = (Particle *)malloc(local_nBodies * sizeof(Particle));
    if (local_particles == NULL) {
        fprintf(stderr, "[Rank %d] Error: Unable to allocate memory for local particles\n", world_rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Allocate memory for all particles (used for MPI_Allgather)
    // Each process needs a copy of ALL particles for force calculation
    Particle *global_particles = (Particle *)malloc(total_nBodies * sizeof(Particle));
    if (global_particles == NULL) {
        fprintf(stderr, "[Rank %d] Error: Unable to allocate memory for global particles\n", world_rank);
        free(local_particles);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Randomize local particles
    randomizeLocalParticles(local_particles, local_nBodies, world_rank);

    // Initial gather of all particles from all processes
    // This fills global_particles on all processes with the initially generated particles
    // MPI_Allgather's recvcount for the receive buffer must be the size of the data from one process.
    // The actual recvcount for MPI_Allgather should be the maximum number of particles any process sends,
    // which is total_nBodies / world_size + (total_nBodies % world_size if last process).
    // However, it's simpler and safer to just use total_nBodies as the recvcount for global_particles
    // and let MPI handle the block distribution based on sendcounts if using MPI_Allgatherv.
    // For MPI_Allgather, the recvcount is the amount of data *each* process receives from *each* process.
    // So, we need to consider the block size.
    // A simpler approach for MPI_Allgather is to ensure uniform block size.
    // The previous code had a slight mismatch in the MPI_Allgather recvcount, it should be the size of the block for *one* process.
    // Corrected to use a consistent block size for MPI_Allgather's recvcount.
    // MPI_Allgather requires that each process sends and receives the same amount of data from each other process.
    // If particles are unevenly distributed (due to total_nBodies % world_size != 0), MPI_Allgather is not the right choice.
    // MPI_Allgatherv would be needed. However, given the context, let's ensure total_nBodies is a multiple of world_size for simplicity
    // or manually handle the last process's contribution during the gather.
    // For this simulation, since computeAccelerations needs all particles, we can stick with MPI_Allgather
    // and let the first `total_nBodies / world_size` elements be filled from each rank, and then handle the remainder.
    // A more robust solution for uneven distribution is MPI_Allgatherv.
    // For now, let's assume `total_nBodies` is divisible by `world_size` or handle the remainder carefully.

    // A correct usage of MPI_Allgather with potentially uneven local_nBodies (due to remainder)
    // would involve each process receiving `total_nBodies / world_size` from all other processes,
    // and the last process receiving `total_nBodies / world_size + (total_nBodies % world_size)` from itself.
    // This is not what MPI_Allgather does. It expects uniform sizes.
    // The most robust way to handle this without `MPI_Allgatherv` (which is more complex) is to pad local_particles
    // or ensure `total_nBodies` is a multiple of `world_size` during testing.

    // For the current structure where `local_nBodies` might vary, `MPI_Allgatherv` is more appropriate.
    // Let's create `sendcounts` and `displs` arrays for `MPI_Allgatherv`.
    int *sendcounts = (int *)malloc(world_size * sizeof(int));
    int *displs = (int *)malloc(world_size * sizeof(int));
    
    // Calculate sendcounts and displacements for MPI_Allgatherv
    for (int i = 0; i < world_size; ++i) {
        sendcounts[i] = total_nBodies / world_size;
        if (i == world_size - 1) {
            sendcounts[i] += (total_nBodies % world_size);
        }
        displs[i] = (i > 0) ? (displs[i-1] + sendcounts[i-1]) : 0;
    }

    // Use MPI_Allgatherv for initial distribution of particles
    MPI_Allgatherv(local_particles, local_nBodies, MPI_Particle_Type,
                   global_particles, sendcounts, displs, MPI_Particle_Type,
                   MPI_COMM_WORLD);

    // Free the temporary sendcounts and displs after initial Allgatherv
    free(sendcounts);
    free(displs);


    // Initialize accelerations to zero for local particles (redundant if randomized to 0, but good for clarity)
    for (int i = 0; i < local_nBodies; i++) {
        local_particles[i].ax = local_particles[i].ay = local_particles[i].az = 0.0f;
        local_particles[i].ax_prev = local_particles[i].ay_prev = local_particles[i].az_prev = 0.0f;
    }
    
    // Compute initial accelerations
    computeAccelerations(local_particles, global_particles, local_nBodies, total_nBodies, start_idx);

    if (world_rank == 0) {
        printf("\n=== N-Body Simulation Results (MPI Velocity Verlet) ===\n");
        printf("Total number of particles: %d\n", total_nBodies);
        printf("Time step: %f\n", dt);
        printf("Number of iterations: %d\n", nIters);
        printf("Integration method: Velocity Verlet (MPI)\n");
        printf("Number of MPI processes: %d\n", world_size);
        printf("=========================================================\n");
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before timing
    start_total = MPI_Wtime();

    for (int iter = 1; iter <= nIters; iter++) {
        start_iter = MPI_Wtime();

        // Perform Verlet integration step
        // This function will also call computeAccelerations internally
        verletIntegration(local_particles, global_particles, dt, local_nBodies, total_nBodies, start_idx);

        end_iter = MPI_Wtime();
        if (world_rank == 0) {
            printf("Iteration %d of %d completed in %f seconds\n", iter, nIters, end_iter - start_iter);
        }
    }

    end_total = MPI_Wtime();
    double total_elapsed_time = end_total - start_total;

    // Collect total simulation time from all processes to find max/min/avg, here just rank 0 gets the total
    double max_total_time;
    MPI_Reduce(&total_elapsed_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Gather all particle states to rank 0 for writing to a single file
    // Need to re-calculate sendcounts and displs for MPI_Gatherv for the final gather
    sendcounts = (int *)malloc(world_size * sizeof(int)); // Reuse variable name, allocate new memory
    displs = (int *)malloc(world_size * sizeof(int));     // Reuse variable name, allocate new memory

    for (int i = 0; i < world_size; ++i) {
        sendcounts[i] = total_nBodies / world_size;
        if (i == world_size - 1) {
            sendcounts[i] += (total_nBodies % world_size);
        }
        displs[i] = (i > 0) ? (displs[i-1] + sendcounts[i-1]) : 0;
    }

    if (world_rank == 0) {
        MPI_Gatherv(local_particles, local_nBodies, MPI_Particle_Type,
                    global_particles, sendcounts, displs, MPI_Particle_Type,
                    0, MPI_COMM_WORLD);

        // Write the output to a file
        FILE *fileWrite = fopen("parallel_output_mpi.txt", "wb");
        if (fileWrite != NULL) {
            size_t written = fwrite(global_particles, sizeof(Particle), total_nBodies, fileWrite);
            fclose(fileWrite);
            if (written == total_nBodies) {
                printf("Final particle states written to parallel_output_mpi.txt\n");
            } else {
                printf("Warning: Could not write all particles to output file\n");
            }
        } else {
            printf("Warning: Could not create output file parallel_output_mpi.txt\n");
        }
    } else {
        // Other processes send their data
        MPI_Gatherv(local_particles, local_nBodies, MPI_Particle_Type,
                    NULL, NULL, NULL, MPI_Particle_Type,
                    0, MPI_COMM_WORLD);
    }
    
    // Free sendcounts and displs after final Gatherv
    free(sendcounts);
    free(displs);


    if (world_rank == 0) {
        double avg_time_per_iter = max_total_time / (double)nIters;
        printf("\n=== Performance Summary ===\n");
        printf("Avg iteration time (based on max total time): %f seconds\n", avg_time_per_iter);
        printf("Total simulation time (max across processes): %f seconds\n", max_total_time);
        printf("Particles processed (total): %d\n", total_nBodies);
        printf("Integration method: Velocity Verlet (MPI)\n");
        printf("Number of MPI processes: %d\n", world_size);
        printf("==========================\n\n");
    }

    // Free custom MPI datatype
    MPI_Type_free(&MPI_Particle_Type);

    // Free allocated memory
    free(local_particles);
    free(global_particles);

    // Finalize MPI environment
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

/* Create custom MPI_Datatype for Particle struct */
void createMPITypeForParticle() {
    const int num_fields = 13; // mass, x,y,z, vx,vy,vz, ax,ay,az, ax_prev,ay_prev,az_prev
    MPI_Aint offsets[num_fields];
    int blocklengths[num_fields];
    MPI_Datatype types[num_fields];

    Particle p_dummy; // Dummy particle to get offsets of members

    // Define types and blocklengths for each member
    for (int i = 0; i < num_fields; ++i) {
        blocklengths[i] = 1;
        types[i] = MPI_FLOAT; // All members are floats
    }

    // Get byte offsets of each member from the beginning of the struct
    MPI_Get_address(&p_dummy.mass, &offsets[0]);
    MPI_Get_address(&p_dummy.x, &offsets[1]);
    MPI_Get_address(&p_dummy.y, &offsets[2]);
    MPI_Get_address(&p_dummy.z, &offsets[3]);
    MPI_Get_address(&p_dummy.vx, &offsets[4]);
    MPI_Get_address(&p_dummy.vy, &offsets[5]);
    MPI_Get_address(&p_dummy.vz, &offsets[6]);
    MPI_Get_address(&p_dummy.ax, &offsets[7]);
    MPI_Get_address(&p_dummy.ay, &offsets[8]);
    MPI_Get_address(&p_dummy.az, &offsets[9]);
    MPI_Get_address(&p_dummy.ax_prev, &offsets[10]);
    MPI_Get_address(&p_dummy.ay_prev, &offsets[11]);
    MPI_Get_address(&p_dummy.az_prev, &offsets[12]);

    // Calculate relative offsets
    MPI_Aint base_address;
    MPI_Get_address(&p_dummy, &base_address);
    for (int i = 0; i < num_fields; ++i) {
        offsets[i] -= base_address;
    }

    // Create the structured MPI datatype
    MPI_Type_create_struct(num_fields, blocklengths, offsets, types, &MPI_Particle_Type);
    MPI_Type_commit(&MPI_Particle_Type);
}

/* Initialize local particle states with random values */
void randomizeLocalParticles(Particle *particles, int n_local, int rank) {
    // Use a unique seed for each process based on current time and rank
    srand(time(NULL) + rank);
    
    for (int i = 0; i < n_local; i++) {
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
}

/* Function that computes accelerations for local particles using all global particles */
void computeAccelerations(Particle *local_p, Particle *global_p, int local_n, int total_n, int start_idx) {
    // Reset accelerations for local particles
    for (int i = 0; i < local_n; i++) {
        local_p[i].ax = local_p[i].ay = local_p[i].az = 0.0f;
    }
    
    // Compute gravitational forces and accelerations
    // Each local particle interacts with ALL global particles
    for (int i = 0; i < local_n; i++) {
        float ax_local = 0.0f, ay_local = 0.0f, az_local = 0.0f;
        
        for (int j = 0; j < total_n; j++) {
            // Check to avoid computing force of a particle on itself
            // local_p[i] corresponds to global_p[start_idx + i]
            if ((start_idx + i) != j) {
                float dx = global_p[j].x - local_p[i].x;
                float dy = global_p[j].y - local_p[i].y;
                float dz = global_p[j].z - local_p[i].z;
                float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
                float invDist = 1.0f / sqrtf(distSqr);
                float invDist3 = invDist * invDist * invDist;

                // F = G * m1 * m2 / r^2, normalize G=1 for simplicity
                // a = F/m = G * m_other / r^2
                float acceleration = global_p[j].mass * invDist3;
                
                ax_local += acceleration * dx;
                ay_local += acceleration * dy;
                az_local += acceleration * dz;
            }
        }
        
        // Store the computed accelerations
        local_p[i].ax = ax_local;
        local_p[i].ay = ay_local;
        local_p[i].az = az_local;
    }
}

/* Velocity Verlet integration step with MPI parallelization */
void verletIntegration(Particle *local_p, Particle *global_p, float dt, int local_n, int total_n, int start_idx) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get world_size correctly inside the function

    // Store current accelerations as previous
    for (int i = 0; i < local_n; i++) {
        local_p[i].ax_prev = local_p[i].ax;
        local_p[i].ay_prev = local_p[i].ay;
        local_p[i].az_prev = local_p[i].az;
    }
    
    // Update positions: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    for (int i = 0; i < local_n; i++) {
        local_p[i].x += local_p[i].vx * dt + 0.5f * local_p[i].ax * dt * dt;
        local_p[i].y += local_p[i].vy * dt + 0.5f * local_p[i].ay * dt * dt;
        local_p[i].z += local_p[i].vz * dt + 0.5f * local_p[i].az * dt * dt;
    }
    
    // All processes need to update their 'global_particles' array with the new positions
    // before computing new accelerations.
    // Re-calculate sendcounts and displs for MPI_Allgatherv, as they might be freed outside this function scope
    int *sendcounts = (int *)malloc(world_size * sizeof(int));
    int *displs = (int *)malloc(world_size * sizeof(int));
    
    for (int i = 0; i < world_size; ++i) {
        sendcounts[i] = total_n / world_size;
        if (i == world_size - 1) {
            sendcounts[i] += (total_n % world_size);
        }
        displs[i] = (i > 0) ? (displs[i-1] + sendcounts[i-1]) : 0;
    }

    MPI_Allgatherv(local_p, local_n, MPI_Particle_Type,
                  global_p, sendcounts, displs, MPI_Particle_Type,
                  MPI_COMM_WORLD);

    free(sendcounts);
    free(displs);

    // Compute new accelerations at t+dt using the updated global positions
    computeAccelerations(local_p, global_p, local_n, total_n, start_idx);
    
    // Update velocities: v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]*dt
    for (int i = 0; i < local_n; i++) {
        local_p[i].vx += 0.5f * (local_p[i].ax_prev + local_p[i].ax) * dt;
        local_p[i].vy += 0.5f * (local_p[i].ay_prev + local_p[i].ay) * dt;
        local_p[i].vz += 0.5f * (local_p[i].az_prev + local_p[i].az) * dt;
    }
}
