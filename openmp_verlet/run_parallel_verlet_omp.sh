#!/bin/bash

# Simple OpenMP N-Body Benchmark Script
PROGRAM_NAME="parallel_nBody_omp"
PARTICLE_GEN="particle_production_omp"
OUTPUT_FILE="parallel_omp_results.txt"
CSV_FILE="parallel_timing_results.csv"

# Thread configurations to test
THREAD_CONFIGS=(4 8 16 32)

echo "=== N-Body OpenMP Benchmark ==="
echo "Date: $(date)"
echo "Thread configs: ${THREAD_CONFIGS[*]}"
echo "================================"

# Clean up existing files
rm -f "$PROGRAM_NAME" "$PARTICLE_GEN" "$OUTPUT_FILE" "$CSV_FILE"
rm -f particles.txt parallel_output.txt

# Compile programs
echo "Compiling..."
gcc -fopenmp -o "$PROGRAM_NAME" parallel_nBody_omp.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile $PROGRAM_NAME"
    exit 1
fi

gcc -fopenmp -o "$PARTICLE_GEN" particle_production_omp.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile $PARTICLE_GEN"
    exit 1
fi

echo "Compilation successful."

# Initialize output files
echo "=== N-Body OpenMP Benchmark Results ===" > "$OUTPUT_FILE"
echo "Date: $(date)" >> "$OUTPUT_FILE"
echo "=======================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Initialize CSV
echo "Particles,Threads,AvgTime,TotalTime" > "$CSV_FILE"

# Function to run simulation
run_simulation() {
    local particles=$1
    local threads=$2
    
    echo "Running: $particles particles, $threads threads..."
    
    # Set OpenMP threads
    export OMP_NUM_THREADS=$threads
    
    # Run simulation and capture output
    output=$(./parallel_nBody_omp $particles 2>&1)
    
    # Extract timing data
    avg_time=$(echo "$output" | grep "Avg iteration time:" | awk '{print $4}')
    total_time=$(echo "$output" | grep "Total simulation time:" | awk '{print $4}')
    
    # Check if we got valid timing data
    if [ -z "$avg_time" ] || [ -z "$total_time" ]; then
        echo "Warning: Could not extract timing data for $particles particles, $threads threads"
        avg_time="ERROR"
        total_time="ERROR"
    fi
    
    # Write to output file
    echo "Particles: $particles, Threads: $threads" >> "$OUTPUT_FILE"
    echo "Avg time: $avg_time sec, Total time: $total_time sec" >> "$OUTPUT_FILE"
    echo "---" >> "$OUTPUT_FILE"
    
    # Write to CSV
    echo "$particles,$threads,$avg_time,$total_time" >> "$CSV_FILE"
    
    # Clean up temp files
    rm -f particles.txt parallel_output.txt
}

# Main benchmark loop
echo "Starting benchmarks..."

# Test particle counts: 1K-10K (steps of 1K), then 20K-100K (steps of 10K)
PARTICLE_COUNTS=()

# Add 1K to 10K
for i in $(seq 1000 1000 10000); do
    PARTICLE_COUNTS+=($i)
done

# Add 20K to 100K
for i in $(seq 20000 10000 100000); do
    PARTICLE_COUNTS+=($i)
done

# Run benchmarks for each thread configuration
for threads in "${THREAD_CONFIGS[@]}"; do
    echo ""
    echo "=== Testing with $threads threads ==="
    
    for particles in "${PARTICLE_COUNTS[@]}"; do
        run_simulation $particles $threads
    done
done

# Summary
echo ""
echo "=== BENCHMARK COMPLETED ==="
echo "Results saved to:"
echo "  - $OUTPUT_FILE (detailed results)"
echo "  - $CSV_FILE (CSV format)"
echo ""
echo "Total simulations run: $((${#PARTICLE_COUNTS[@]} * ${#THREAD_CONFIGS[@]}))"
echo "Particle counts tested: ${PARTICLE_COUNTS[*]}"
echo "Thread configurations: ${THREAD_CONFIGS[*]}"

# Clean up executables
rm -f "$PROGRAM_NAME" "$PARTICLE_GEN"

echo "Done!"