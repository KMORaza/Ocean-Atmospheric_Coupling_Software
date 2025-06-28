import numpy as np

class VariableResolutionGrid:
    def __init__(self, grid_size, coast_factor):
        self.grid_size = grid_size
        self.coast_factor = coast_factor  # Factor for finer resolution near coasts
        
    def get_spatial_steps(self):
        # Simulate finer resolution near "coasts" (edges) and coarser in the center
        dx = np.ones((self.grid_size, self.grid_size))
        dy = np.ones((self.grid_size, self.grid_size))
        
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                # Finer resolution near edges (simulating coasts)
                if i < 5 or i >= self.grid_size-5 or j < 5 or j >= self.grid_size-5:
                    dx[i,j] = 1.0 / self.coast_factor
                    dy[i,j] = 1.0 / self.coast_factor
                else:
                    dx[i,j] = 1.0
                    dy[i,j] = 1.0
        return dx, dy