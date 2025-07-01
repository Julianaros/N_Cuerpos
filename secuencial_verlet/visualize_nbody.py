import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation, PillowWriter
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')

class NBodyVisualizer:
    def __init__(self, positions_file='particle_positions.csv', velocity_file='velocity_field.csv'):
        """
        Initialize the visualizer with data files
        """
        self.positions_file = positions_file
        self.velocity_file = velocity_file
        self.positions_data = None
        self.velocity_data = None
        
    def load_data(self):
        """
        Load the simulation data from CSV files
        """
        try:
            print("Loading particle position data...")
            self.positions_data = pd.read_csv(self.positions_file)
            print(f"Loaded {len(self.positions_data)} position records")
            
            print("Loading velocity field data...")
            self.velocity_data = pd.read_csv(self.velocity_file)
            print(f"Loaded {len(self.velocity_data)} velocity records")
            
            print("Data loading complete!")
            
        except FileNotFoundError as e:
            print(f"Error: Could not find data file - {e}")
            print("Make sure to run the modified C simulation first!")
            return False
        
        return True
    
    def create_density_heatmap_gif(self, projection='xy', bins=50, start_iter=0, end_iter=50, 
                                   interval=200, save_gif=True, filename=None):
        """
        Create an animated GIF of density heatmap evolution
        
        Parameters:
        - projection: 'xy', 'xz', or 'yz' plane
        - bins: number of bins for the histogram
        - start_iter: starting iteration
        - end_iter: ending iteration
        - interval: milliseconds between frames
        - save_gif: whether to save the GIF
        - filename: custom filename for the GIF
        """
        if self.positions_data is None:
            print("Error: No position data loaded!")
            return
        
        # Select projection coordinates
        if projection == 'xy':
            x_col, y_col = 'x', 'y'
            xlabel, ylabel = 'X Position', 'Y Position'
        elif projection == 'xz':
            x_col, y_col = 'x', 'z'
            xlabel, ylabel = 'X Position', 'Z Position'
        elif projection == 'yz':
            x_col, y_col = 'y', 'z'
            xlabel, ylabel = 'Y Position', 'Z Position'
        else:
            print("Error: projection must be 'xy', 'xz', or 'yz'")
            return
        
        # Get available iterations in the specified range
        available_iters = sorted(self.positions_data['iteration'].unique())
        iterations_to_plot = [iter for iter in available_iters if start_iter <= iter <= end_iter]
        
        if not iterations_to_plot:
            print(f"No iterations found in range {start_iter}-{end_iter}")
            return
        
        print(f"Creating GIF with {len(iterations_to_plot)} frames: iterations {iterations_to_plot[0]} to {iterations_to_plot[-1]}")
        
        # Use fixed bounds for consistent scaling (-25 to 25)
        x_min, x_max = -25.0, 25.0
        y_min, y_max = -25.0, 25.0
        
        # Calculate global max density for consistent color scaling
        max_density = 0
        for iteration in iterations_to_plot:
            iter_data = self.positions_data[self.positions_data['iteration'] == iteration]
            if len(iter_data) > 0:
                x_data = iter_data[x_col].values
                y_data = iter_data[y_col].values
                hist, _, _ = np.histogram2d(x_data, y_data, bins=bins, 
                                         range=[[x_min, x_max], [y_min, y_max]])
                max_density = max(max_density, hist.max())
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(10, 8))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')
        
        # Initialize the image with fixed bounds (-25 to 25)
        x_min, x_max = -25.0, 25.0
        y_min, y_max = -25.0, 25.0
        
        im = ax.imshow(np.zeros((bins, bins)), origin='lower', cmap='hot', 
                      extent=[x_min, x_max, y_min, y_max], aspect='auto',
                      vmin=0, vmax=max_density, interpolation='bilinear')
        
        # Set up the plot
        ax.set_xlabel(xlabel, fontsize=12, color='white')
        ax.set_ylabel(ylabel, fontsize=12, color='white')
        ax.tick_params(colors='white')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Particle Density', fontsize=12, color='white')
        cbar.ax.tick_params(colors='white')
        
        # Title and info text
        title = ax.set_title(f'Particle Density Evolution - {projection.upper()} Plane\nIteration: 0', 
                           fontsize=14, fontweight='bold', color='white', pad=20)
        
        # Info text box
        info_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, fontsize=10,
                           verticalalignment='top', color='white',
                           bbox=dict(boxstyle='round', facecolor='black', alpha=0.7, edgecolor='white'))
        
        plt.tight_layout()
        
        def animate(frame_idx):
            """Animation function"""
            iteration = iterations_to_plot[frame_idx]
            
            # Get data for this iteration
            iter_data = self.positions_data[self.positions_data['iteration'] == iteration]
            
            if len(iter_data) == 0:
                return [im, title, info_text]
            
            # Create histogram
            x_data = iter_data[x_col].values
            y_data = iter_data[y_col].values
            hist, _, _ = np.histogram2d(x_data, y_data, bins=bins, 
                                     range=[[x_min, x_max], [y_min, y_max]])
            
            # Update image
            im.set_array(hist.T)
            
            # Update title
            title.set_text(f'Particle Density Evolution - {projection.upper()} Plane\nIteration: {iteration}')
            
            # Update info text
            particles_count = len(iter_data)
            max_local_density = hist.max()
            progress = (frame_idx + 1) / len(iterations_to_plot) * 100
            
            info_str = f'Particles: {particles_count}\nMax Density: {max_local_density:.0f}\nProgress: {progress:.1f}%'
            info_text.set_text(info_str)
            
            return [im, title, info_text]
        
        # Create animation
        print("Creating animation frames...")
        anim = FuncAnimation(fig, animate, frames=len(iterations_to_plot), 
                           interval=interval, blit=False, repeat=True)
        
        # Save GIF
        if save_gif:
            if filename is None:
                filename = f'density_heatmap_evolution_{projection}_{start_iter}to{end_iter}.gif'
            
            print(f"Saving GIF as {filename}...")
            print("This may take a while depending on the number of frames...")
            
            # Use PillowWriter for better GIF quality
            writer = PillowWriter(fps=1000//interval)  # Convert interval to fps
            anim.save(filename, writer=writer, dpi=150)
            
            print(f"GIF saved successfully as {filename}")
        
        plt.show()
        return anim
    
    def create_density_heatmap(self, iteration=0, projection='xy', bins=50, save_fig=True):
        """
        Create a 2D density heatmap of particle positions
        
        Parameters:
        - iteration: which simulation iteration to plot
        - projection: 'xy', 'xz', or 'yz' plane
        - bins: number of bins for the histogram
        - save_fig: whether to save the figure
        """
        if self.positions_data is None:
            print("Error: No position data loaded!")
            return
        
        # Filter data for specific iteration
        iter_data = self.positions_data[self.positions_data['iteration'] == iteration]
        
        if len(iter_data) == 0:
            print(f"Warning: No data found for iteration {iteration}")
            return
        
        # Select projection coordinates
        if projection == 'xy':
            x_col, y_col = 'x', 'y'
            xlabel, ylabel = 'X Position', 'Y Position'
        elif projection == 'xz':
            x_col, y_col = 'x', 'z'
            xlabel, ylabel = 'X Position', 'Z Position'
        elif projection == 'yz':
            x_col, y_col = 'y', 'z'
            xlabel, ylabel = 'Y Position', 'Z Position'
        else:
            print("Error: projection must be 'xy', 'xz', or 'yz'")
            return
        
        # Create the heatmap
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create 2D histogram
        x_data = iter_data[x_col].values
        y_data = iter_data[y_col].values
        
        # Use fixed bounds (-25 to 25) for consistent scaling
        x_min, x_max = -25.0, 25.0
        y_min, y_max = -25.0, 25.0
        
        # Create histogram
        hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=bins, 
                                                range=[[x_min, x_max], [y_min, y_max]])
        
        # Plot heatmap
        im = ax.imshow(hist.T, origin='lower', cmap='hot', interpolation='bilinear',
                      extent=[x_min, x_max, y_min, y_max], aspect='auto')
        
        # Customize plot
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(f'Particle Density Heatmap - {projection.upper()} Plane\nIteration {iteration}', 
                    fontsize=14, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Particle Count', fontsize=12)
        
        # Add statistics text
        stats_text = f'Total Particles: {len(iter_data)}\nBins: {bins}×{bins}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        if save_fig:
            filename = f'density_heatmap_{projection}_iter_{iteration:03d}.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Density heatmap saved as {filename}")
        
        plt.show()
    
    def create_velocity_magnitude_map(self, iteration=0, projection='xy', save_fig=True):
        """
        Create a 2D velocity magnitude map (solo magnitud, sin vectores)
        
        Parameters:
        - iteration: which simulation iteration to plot
        - projection: 'xy', 'xz', or 'yz' plane
        - save_fig: whether to save the figure
        """
        if self.velocity_data is None:
            print("Error: No velocity data loaded!")
            return
        
        # Filter data for specific iteration
        iter_data = self.velocity_data[self.velocity_data['iteration'] == iteration]
        
        if len(iter_data) == 0:
            print(f"Warning: No data found for iteration {iteration}")
            return
        
        # Select projection coordinates
        if projection == 'xy':
            x_col, y_col = 'x', 'y'
            xlabel, ylabel = 'X Position', 'Y Position'
        elif projection == 'xz':
            x_col, y_col = 'x', 'z'
            xlabel, ylabel = 'X Position', 'Z Position'
        elif projection == 'yz':
            x_col, y_col = 'y', 'z'
            xlabel, ylabel = 'Y Position', 'Z Position'
        else:
            print("Error: projection must be 'xy', 'xz', or 'yz'")
            return
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Get data
        x_data = iter_data[x_col].values
        y_data = iter_data[y_col].values
        v_mag = iter_data['velocity_magnitude'].values
        
        # Create scatter plot with velocity magnitude as color
        scatter = ax.scatter(x_data, y_data, c=v_mag, cmap='viridis', 
                             s=30, alpha=0.8, edgecolors='none')
        
        # Set fixed axis limits (-25 to 25)
        ax.set_xlim(-25, 25)
        ax.set_ylim(-25, 25)
        
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(f'Velocity Magnitude Field - {projection.upper()} Plane\nIteration {iteration}', 
                     fontsize=14, fontweight='bold')
        
        # Add colorbar for magnitude
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Velocity Magnitude', fontsize=12)
        
        # Add statistics
        stats_text = f'Particles: {len(iter_data)}\nMax vel: {v_mag.max():.3f}\nAvg vel: {v_mag.mean():.3f}\nMin vel: {v_mag.min():.3f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        if save_fig:
            filename = f'velocity_magnitude_{projection}_iter_{iteration:03d}.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Velocity magnitude map saved as {filename}")
        
        plt.show()
    
    def create_evolution_animation(self, projection='xy', max_frames=10):
        """
        Create evolution animation - placeholder implementation
        """
        print("Evolution animation functionality not implemented in this version.")
        pass
    
    def generate_summary_report(self):
        """
        Generate a summary report of the simulation data
        """
        if self.positions_data is None or self.velocity_data is None:
            print("Error: Data not loaded!")
            return
        
        print("\n" + "="*60)
        print("N-BODY SIMULATION ANALYSIS REPORT")
        print("="*60)
        
        # Basic statistics
        n_particles = self.positions_data['particle_id'].nunique()
        n_iterations = self.positions_data['iteration'].nunique()
        
        print(f"Number of particles: {n_particles}")
        print(f"Number of iterations: {n_iterations}")
        print(f"Total data points: {len(self.positions_data)}")
        
        # Spatial bounds
        print(f"\nSpatial bounds:")
        for coord in ['x', 'y', 'z']:
            min_val = self.positions_data[coord].min()
            max_val = self.positions_data[coord].max()
            print(f"  {coord.upper()}: [{min_val:.3f}, {max_val:.3f}]")
        
        # Velocity statistics
        print(f"\nVelocity statistics:")
        v_mag = self.velocity_data['velocity_magnitude']
        print(f"  Maximum velocity: {v_mag.max():.6f}")
        print(f"  Average velocity: {v_mag.mean():.6f}")
        print(f"  Velocity std dev: {v_mag.std():.6f}")
        
        # Final positions vs initial
        initial_data = self.positions_data[self.positions_data['iteration'] == 0]
        final_data = self.positions_data[self.positions_data['iteration'] == self.positions_data['iteration'].max()]
        
        initial_center = np.array([initial_data['x'].mean(), initial_data['y'].mean(), initial_data['z'].mean()])
        final_center = np.array([final_data['x'].mean(), final_data['y'].mean(), final_data['z'].mean()])
        center_movement = np.linalg.norm(final_center - initial_center)
        
        print(f"\nSystem evolution:")
        print(f"  Initial center of mass: ({initial_center[0]:.3f}, {initial_center[1]:.3f}, {initial_center[2]:.3f})")
        print(f"  Final center of mass: ({final_center[0]:.3f}, {final_center[1]:.3f}, {final_center[2]:.3f})")
        print(f"  Center of mass movement: {center_movement:.6f}")
        
        print("="*60)

def main():
    """
    Main function to demonstrate the visualization capabilities
    """
    print("N-Body Simulation Visualizer with GIF Animation")
    print("="*50)
    
    # Initialize visualizer
    viz = NBodyVisualizer()
    
    # Load data
    if not viz.load_data():
        return
    
    # Generate report
    viz.generate_summary_report()
    
    print("\nGenerating visualizations...")
    
    # Create the main GIF animation of density heatmap evolution
    print("\n*** Creating GIF Animation of Density Heatmap Evolution ***")
    print("This will create a GIF showing particle concentration from iteration 0 to 50")
    
    # Create GIF for XY plane (most common view)
    viz.create_density_heatmap_gif(
        projection='xy',
        bins=60,  # Higher resolution for better detail
        start_iter=0,
        end_iter=50,
        interval=300,  # 300ms between frames (slower for better viewing)
        save_gif=True,
        filename='particle_density_evolution_0to50.gif'
    )
    
    # Optional: Create GIFs for other projections
    create_other_projections = input("\nDo you want to create GIFs for XZ and YZ projections too? (y/n): ").lower().strip()
    
    if create_other_projections == 'y':
        print("Creating XZ projection GIF...")
        viz.create_density_heatmap_gif(
            projection='xz',
            bins=60,
            start_iter=0,
            end_iter=50,
            interval=300,
            save_gif=True,
            filename='particle_density_evolution_xz_0to50.gif'
        )
        
        print("Creating YZ projection GIF...")
        viz.create_density_heatmap_gif(
            projection='yz',
            bins=60,
            start_iter=0,
            end_iter=50,
            interval=300,
            save_gif=True,
            filename='particle_density_evolution_yz_0to50.gif'
        )
    
    # Create some static heatmaps for comparison
    print("\nCreating static density heatmaps for key iterations...")
    key_iterations = [0, 10, 25, 50]
    
    for iteration in key_iterations:
        if iteration in viz.positions_data['iteration'].values:
            print(f"Creating static heatmap for iteration {iteration}")
            viz.create_density_heatmap(iteration=iteration, projection='xy')
    
    print("\n" + "="*50)
    print("VISUALIZATION COMPLETE!")
    print("="*50)
    print("Generated files:")
    print("- particle_density_evolution_0to50.gif (main GIF animation)")
    if create_other_projections == 'y':
        print("- particle_density_evolution_xz_0to50.gif")
        print("- particle_density_evolution_yz_0to50.gif")
    print("- Various static PNG heatmaps")
    print("\nThe GIF shows the evolution of particle concentration over time.")
    print("Hot colors (red/yellow) indicate high particle density areas.")

if __name__ == "__main__":
    main()