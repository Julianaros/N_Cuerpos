#!/bin/bash
# MPI N-Body Benchmark Script - Fixed Version with OpenMPI Oversubscribe
PROGRAM_NAME="parallel_nBody_mpi"
OUTPUT_FILE="parallel_mpi_results.txt"
CSV_FILE="parallel_timing_results.csv"

# Number of MPI processes to test
PROCESS_CONFIGS=(4 8 16 32) # You can adjust these based on your system's core count

echo "=== N-Body MPI Benchmark ==="
echo "Date: $(date)"
echo "Process configs: ${PROCESS_CONFIGS[*]}"
echo "================================"

# Clean up existing files
rm -f "$PROGRAM_NAME" "$OUTPUT_FILE" "$CSV_FILE"
rm -f parallel_output_mpi.txt

# Compile program
echo "Compiling..."
mpicc -o "$PROGRAM_NAME" parallel_nBody_mpi.c -lm -O2
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile $PROGRAM_NAME"
    exit 1
fi
echo "Compilation successful."

# Initialize output files
echo "=== N-Body MPI Benchmark Results ===" > "$OUTPUT_FILE"
echo "Date: $(date)" >> "$OUTPUT_FILE"
echo "=======================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Initialize CSV
echo "Particles,Processes,AvgTime,TotalTime" > "$CSV_FILE"

# Function to run simulation
run_simulation() {
    local particles=$1
    local processes=$2
    
    echo "Running: $particles particles, $processes processes..."
    
    # Determine MPI options based on process count
    local mpi_options=""
    if [ "$processes" -gt 16 ]; then
        # For more than 16 processes, use oversubscribe to allow hyperthreading
        mpi_options="--oversubscribe"
        echo "  Using --oversubscribe for $processes processes (more than physical cores)"
    fi
    
    # Run simulation and capture output from rank 0
    # Use timeout to prevent hanging and redirect stderr to stdout
    output=$(timeout 300s mpirun -np "$processes" $mpi_options "./$PROGRAM_NAME" "$particles" 2>&1)
    exit_code=$?
    
    # Check if simulation completed successfully
    if [ $exit_code -ne 0 ]; then
        echo "Warning: Simulation failed or timed out for $particles particles, $processes processes"
        echo "Exit code: $exit_code"
        # Show some debug info for failed runs
        if [ $exit_code -eq 124 ]; then
            echo "  -> Reason: Timeout (>300s)"
        elif [ $exit_code -eq 1 ]; then
            echo "  -> Reason: Runtime error or MPI issue"
        fi
        avg_time="TIMEOUT/ERROR"
        total_time="TIMEOUT/ERROR"
    else
        # Extract timing data from the output using sed (more reliable than awk for this)
        # Expected format:
        # "Avg iteration time (based on max total time): X.XXXXXX seconds"
        # "Total simulation time (max across processes): X.XXXXXX seconds"
        
        avg_time=$(echo "$output" | sed -n 's/.*Avg iteration time (based on max total time): \([0-9][0-9]*\.[0-9]*\) seconds.*/\1/p')
        total_time=$(echo "$output" | sed -n 's/.*Total simulation time (max across processes): \([0-9][0-9]*\.[0-9]*\) seconds.*/\1/p')
        
        # If sed didn't work, try alternative patterns (in case format is slightly different)
        if [ -z "$avg_time" ]; then
            # Try extracting any decimal number after the colon
            avg_time=$(echo "$output" | grep "Avg iteration time" | sed 's/.*: \([0-9][0-9]*\.[0-9]*\).*/\1/')
        fi
        
        if [ -z "$total_time" ]; then
            # Try extracting any decimal number after the colon
            total_time=$(echo "$output" | grep "Total simulation time" | sed 's/.*: \([0-9][0-9]*\.[0-9]*\).*/\1/')
        fi
        
        # Check if we got valid timing data
        if [ -z "$avg_time" ] || [ -z "$total_time" ]; then
            echo "Warning: Could not extract timing data for $particles particles, $processes processes"
            echo "Debug - Lines containing 'time':"
            echo "$output" | grep -i "time" | head -5
            echo "Debug - Performance section:"
            echo "$output" | grep -A 10 "Performance Summary"
            avg_time="PARSE_ERROR"
            total_time="PARSE_ERROR"
        else
            # Validate that extracted values are numeric
            if ! [[ "$avg_time" =~ ^[0-9]+\.?[0-9]*$ ]]; then
                echo "Warning: Invalid avg_time value: '$avg_time'"
                echo "Debug - Avg time line:"
                echo "$output" | grep "Avg iteration time"
                avg_time="INVALID"
            fi
            if ! [[ "$total_time" =~ ^[0-9]+\.?[0-9]*$ ]]; then
                echo "Warning: Invalid total_time value: '$total_time'"
                echo "Debug - Total time line:"
                echo "$output" | grep "Total simulation time"
                total_time="INVALID"
            fi
        fi
    fi
    
    # Display results
    echo "  -> Avg time: $avg_time sec, Total time: $total_time sec"
    
    # Write to output file
    echo "Particles: $particles, Processes: $processes" >> "$OUTPUT_FILE"
    echo "Avg time: $avg_time sec, Total time: $total_time sec" >> "$OUTPUT_FILE"
    if [ "$processes" -gt 16 ]; then
        echo "Note: Used --oversubscribe (hyperthreading)" >> "$OUTPUT_FILE"
    fi
    echo "---" >> "$OUTPUT_FILE"
    
    # Write to CSV
    echo "$particles,$processes,$avg_time,$total_time" >> "$CSV_FILE"
    
    # Clean up temp files if created
    rm -f parallel_output_mpi.txt
    
    # Small delay to prevent overwhelming the system
    sleep 1
}

# Main benchmark loop
echo "Starting benchmarks..."

# Test particle counts: 1K-10K (steps of 1K), then 10K-100K (steps of 10K)
PARTICLE_COUNTS=()
# Add 1K to 10K in steps of 1K
for i in $(seq 1000 1000 10000); do
    PARTICLE_COUNTS+=($i)
done
# Add 20K to 100K in steps of 10K
for i in $(seq 20000 10000 100000); do
    PARTICLE_COUNTS+=($i)
done

echo "Particle counts to test: ${PARTICLE_COUNTS[*]}"
echo ""
echo "Note: Processes >16 will use --oversubscribe to utilize hyperthreading"
echo ""

# Run benchmarks for each process configuration
for processes in "${PROCESS_CONFIGS[@]}"; do
    echo ""
    echo "=== Testing with $processes processes ==="
    
    for particles in "${PARTICLE_COUNTS[@]}"; do
        run_simulation $particles $processes
    done
done

# Generate summary statistics
echo ""
echo "=== GENERATING SUMMARY ==="

# Create a summary section in the output file
echo "" >> "$OUTPUT_FILE"
echo "=== SUMMARY STATISTICS ===" >> "$OUTPUT_FILE"
echo "Date: $(date)" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Process CSV to generate summary (excluding error entries)
if command -v awk >/dev/null 2>&1; then
    echo "Performance summary by process count:" >> "$OUTPUT_FILE"
    for processes in "${PROCESS_CONFIGS[@]}"; do
        echo "Processes: $processes" >> "$OUTPUT_FILE"
        awk -F',' -v proc="$processes" '
            NR>1 && $2==proc && $3!~/ERROR|TIMEOUT|INVALID|PARSE_ERROR/ {
                sum+=$3; count++; 
                if(min=="" || $3<min) min=$3; 
                if(max=="" || $3>max) max=$3
            } 
            END {
                if(count>0) {
                    printf "  Avg iteration time - Min: %.6f, Max: %.6f, Mean: %.6f (%d samples)\n", min, max, sum/count, count
                } else {
                    printf "  No valid samples found\n"
                }
            }' "$CSV_FILE" >> "$OUTPUT_FILE"
    done
    
    # Add scaling analysis
    echo "" >> "$OUTPUT_FILE"
    echo "=== SCALING ANALYSIS ===" >> "$OUTPUT_FILE"
    echo "Speedup analysis (compared to 4 processes baseline):" >> "$OUTPUT_FILE"
    
    # Get baseline time for 10000 particles with 4 processes
    baseline=$(awk -F',' 'NR>1 && $1==10000 && $2==4 && $3!~/ERROR|TIMEOUT|INVALID|PARSE_ERROR/ {print $3}' "$CSV_FILE")
    
    if [ -n "$baseline" ]; then
        echo "Baseline (10000 particles, 4 processes): $baseline seconds" >> "$OUTPUT_FILE"
        for processes in "${PROCESS_CONFIGS[@]}"; do
            if [ "$processes" -ne 4 ]; then
                current_time=$(awk -F',' -v proc="$processes" 'NR>1 && $1==10000 && $2==proc && $3!~/ERROR|TIMEOUT|INVALID|PARSE_ERROR/ {print $3}' "$CSV_FILE")
                if [ -n "$current_time" ]; then
                    speedup=$(awk "BEGIN {printf \"%.2f\", $baseline/$current_time}")
                    efficiency=$(awk -v s="$speedup" -v p="$processes" "BEGIN {printf \"%.1f\", (s/p)*100}")
                    echo "  $processes processes: ${current_time}s, Speedup: ${speedup}x, Efficiency: ${efficiency}%" >> "$OUTPUT_FILE"
                fi
            fi
        done
    fi
fi

# Summary
echo ""
echo "=== BENCHMARK COMPLETED ==="
echo "Results saved to:"
echo "  - $OUTPUT_FILE (detailed results)"
echo "  - $CSV_FILE (CSV format)"
echo ""
echo "Total simulations run: $((${#PARTICLE_COUNTS[@]} * ${#PROCESS_CONFIGS[@]}))"
echo "Particle counts tested: ${PARTICLE_COUNTS[*]}"
echo "Process configurations: ${PROCESS_CONFIGS[*]}"

# Show a preview of results
echo ""
echo "=== RESULTS PREVIEW ==="
if [ -f "$CSV_FILE" ]; then
    echo "First few lines of CSV results:"
    head -6 "$CSV_FILE"
    echo ""
    echo "Last few lines of CSV results:"
    tail -5 "$CSV_FILE"
fi

# Performance summary
echo ""
echo "=== QUICK PERFORMANCE SUMMARY ==="
if [ -f "$CSV_FILE" ]; then
    echo "Best times for 10000 particles:"
    awk -F',' 'NR>1 && $1==10000 && $3!~/ERROR|TIMEOUT|INVALID|PARSE_ERROR/ {print "  " $2 " processes: " $3 " seconds"}' "$CSV_FILE" | sort -k3 -n
fi

# Clean up executable
rm -f "$PROGRAM_NAME"
echo ""
echo "Done! Check the output files for detailed results."