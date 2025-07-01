#!/bin/bash

echo "=== N-Body Simulation with Visualization Pipeline ==="
echo "Date: $(date)"
echo "======================================================"

# Clean up any existing files
echo "Cleaning up previous files..."
rm -f sequential_verlet particle_production
rm -f particle_positions.csv velocity_field.csv
rm -f particles.txt sequential_verlet_output.txt
rm -f *.png *.gif

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

# Function to run simulation with specified parameters
run_simulation() {
    local nParticles=$1
    echo "Running N-Body simulation with $nParticles particles..."
    echo "This will generate CSV files for visualization..."
    
    # Run the simulation
    ./sequential_verlet "$nParticles"
    
    if [ $? -eq 0 ]; then
        echo "Simulation completed successfully!"
        echo "Generated files:"
        ls -la *.csv 2>/dev/null || echo "No CSV files found"
    else
        echo "Error: Simulation failed"
        return 1
    fi
}

# Check if user provided number of particles as argument
if [ $# -eq 1 ]; then
    nParticles=$1
else
    # Default number of particles for visualization
    # Using smaller number for better visualization performance
    nParticles=2000
fi

echo "Starting simulation with $nParticles particles..."
run_simulation $nParticles

if [ $? -eq 0 ]; then
    echo ""
    echo "=== Running Python Visualization ==="
    
    # Check if Python is available
    if command -v python3 &> /dev/null; then
        echo "Running visualization script..."
        python3 visualize_nbody.py
    elif command -v python &> /dev/null; then
        echo "Running visualization script..."
        python visualize_nbody.py
    else
        echo "Warning: Python not found. Please install Python to run visualizations."
        echo "You can manually run: python3 visualize_nbody.py"
    fi
    
    echo ""
    echo "=== Pipeline Complete ==="
    echo "Generated files:"
    echo "Data files:"
    ls -la *.csv 2>/dev/null
    echo ""
    echo "Visualization files:"
    ls -la *.png *.gif 2>/dev/null
    echo ""
    echo "To re-run only the visualization:"
    echo "  python3 visualize_nbody.py"
    echo ""
    echo "To run with different number of particles:"
    echo "  ./run_visualization.sh <number_of_particles>"
fi