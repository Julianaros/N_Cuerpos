#!/bin/bash

OUTPUT_FILE="sequential_results.txt"
PROGRAM_NAME="sequential"

echo "=== N-Body Sequential Simulation Benchmark ===" > "$OUTPUT_FILE"
echo "Date: $(date)" >> "$OUTPUT_FILE"
echo "Machine: $(uname -a)" >> "$OUTPUT_FILE"
echo "=========================================" >> "$OUTPUT_FILE"
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
gcc -o sequential sequential_nBody.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile sequential_nBody.c"
    exit 1
fi

gcc -o particle_production particle_production.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile particle_production.c"
    exit 1
fi

echo "Compilation successful."
echo ""

# Function to run simulation and extract timing data
run_simulation() {
    local nParticles=$1
    echo "Running simulation with $nParticles particles..."
    
    # Run the program and capture output
    output=$(./sequential "$nParticles" 2>&1)
    
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
    if [ ! -f "timing_summary.csv" ]; then
        echo "Particles,AvgIterationTime,TotalTime" > "timing_summary.csv"
    fi
    echo "$nParticles,$avg_time,$total_time" >> "timing_summary.csv"
}

# Run simulations from 1000 to 10000 particles (increments of 1000)
echo "Starting benchmark runs..."
for nParticles in $(seq 1000 1000 10000); do
    run_simulation $nParticles
    
    # Clean up intermediate files to save space
    if [ -f "particles.txt" ]; then
        rm "particles.txt"
    fi
    if [ -f "sequential_output.txt" ]; then
        rm "sequential_output.txt"
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
    if [ -f "sequential_output.txt" ]; then
        rm "sequential_output.txt"
    fi
done

echo "=== Benchmark Complete ===" >> "$OUTPUT_FILE"
echo "Results saved to: $OUTPUT_FILE" >> "$OUTPUT_FILE"
echo "CSV summary saved to: timing_summary.csv" >> "$OUTPUT_FILE"

echo ""
echo "Sequential benchmark completed successfully!"
echo "Results saved to: $OUTPUT_FILE"
echo "CSV timing summary saved to: timing_summary.csv"
echo ""
echo "To view results:"
echo "  cat $OUTPUT_FILE"
echo "  cat timing_summary.csv"