#!/bin/bash

OUTPUT_FILE="sequential_verlet_results.txt"
PROGRAM_NAME="sequential_verlet"

echo "=== N-Body Sequential Simulation Benchmark (Velocity Verlet) ===" > "$OUTPUT_FILE"
echo "Date: $(date)" >> "$OUTPUT_FILE"
echo "Machine: $(uname -a)" >> "$OUTPUT_FILE"
echo "Integration Method: Velocity Verlet" >> "$OUTPUT_FILE"
echo "=================================================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Clean up any existing executables
if [ -e "$PROGRAM_NAME" ]; then
    rm "$PROGRAM_NAME"
fi
if [ -e "particle_production" ]; then
    rm "particle_production"
fi

# Compile the programs
echo "Compiling programs..."
gcc -o sequential_verlet sequential_nBody_verlet.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile sequential_nBody_verlet.c"
    exit 1
fi

gcc -o particle_production particle_production_verlet.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile particle_production_verlet.c"
    exit 1
fi

echo "Compilation successful."
echo ""

# Function to run simulation and extract timing data
run_simulation() {
    local nParticles=$1
    echo "Running Velocity Verlet simulation with $nParticles particles..."
    
    # Run the program and capture output
    output=$(./sequential_verlet "$nParticles" 2>&1)
    
    # Extract timing information
    avg_time=$(echo "$output" | grep "Avg iteration time:" | awk '{print $4}')
    total_time=$(echo "$output" | grep "Total simulation time:" | awk '{print $4}')
    
    # Write results in a structured format
    echo "Particles: $nParticles" >> "$OUTPUT_FILE"
    echo "Average iteration time: $avg_time seconds" >> "$OUTPUT_FILE"
    echo "Total simulation time: $total_time seconds" >> "$OUTPUT_FILE"
    echo "Full output:" >> "$OUTPUT_FILE"
    echo "$output" >> "$OUTPUT_FILE"
    echo "=========================================" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    
    # Also create a CSV-friendly summary
    if [ ! -f "timing_summary_verlet.csv" ]; then
        echo "Particles,AvgIterationTime,TotalTime,Method" > "timing_summary_verlet.csv"
    fi
    echo "$nParticles,$avg_time,$total_time,Verlet" >> "timing_summary_verlet.csv"
}

# Run simulations from 1000 to 10000 particles (increments of 1000)
echo "Starting Velocity Verlet benchmark runs..."
for nParticles in $(seq 1000 1000 10000); do
    run_simulation $nParticles
    
    # Clean up intermediate files to save space
    if [ -f "particles.txt" ]; then
        rm "particles.txt"
    fi
    if [ -f "sequential_verlet_output.txt" ]; then
        rm "sequential_verlet_output.txt"
    fi
done

# For larger numbers, use bigger increments to avoid extremely long runs
echo "Running larger simulations with 10K increments..."
for nParticles in $(seq 20000 10000 100000); do
    run_simulation $nParticles
    
    # Clean up intermediate files to save space
    if [ -f "particles.txt" ]; then
        rm "particles.txt"
    fi
    if [ -f "sequential_verlet_output.txt" ]; then
        rm "sequential_verlet_output.txt"
    fi
done

echo "=== Velocity Verlet Benchmark Complete ===" >> "$OUTPUT_FILE"
echo "Results saved to: $OUTPUT_FILE" >> "$OUTPUT_FILE"
echo "CSV summary saved to: timing_summary_verlet.csv" >> "$OUTPUT_FILE"

echo ""
echo "Sequential Velocity Verlet benchmark completed successfully!"
echo "Results saved to: $OUTPUT_FILE"
echo "CSV timing summary saved to: timing_summary_verlet.csv"
echo ""
echo "To view results:"
echo "  cat $OUTPUT_FILE"
echo "  cat timing_summary_verlet.csv"