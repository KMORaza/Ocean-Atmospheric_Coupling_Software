import numpy as np

class VariableResolutionGrid:
    def __init__(self, grid_size, coast_factor):
        self.grid_size = grid_size
        self.coast_factor = coast_factor
        
    def get_spatial_steps(self):
        # Create a smooth transition from fine to coarse resolution
        dx = np.ones((self.grid_size, self.grid_size))
        dy = np.ones((self.grid_size, self.grid_size))
        
        # Simulate coastline with a Gaussian decay
        x, y = np.meshgrid(np.arange(self.grid_size), np.arange(self.grid_size))
        center = self.grid_size / 2
        dist = np.sqrt((x - center)**2 + (y - center)**2)
        max_dist = np.sqrt(2) * center
        
        # Finer near coasts (edges), coarser in center
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                factor = 1 + (self.coast_factor - 1) * np.exp(-((dist[i,j]/max_dist)**2))
                dx[i,j] = 1.0 / factor
                dy[i,j] = 1.0 / factor
        return dx, dy