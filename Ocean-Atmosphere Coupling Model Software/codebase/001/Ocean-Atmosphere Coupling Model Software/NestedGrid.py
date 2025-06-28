import numpy as np

class NestedGrid:
    def __init__(self, global_grid_size, nested_grid_size, coast_factor):
        self.global_grid_size = global_grid_size
        self.nested_grid_size = nested_grid_size
        self.coast_factor = coast_factor
        
        # Define nested grid region (center of global grid)
        self.offset = (global_grid_size - nested_grid_size) // 2
        self.nested_dx = 1.0 / (2 * coast_factor)  # Higher resolution
        
    def update(self, To, Ta, k, Co, Ca, dt):
        To_n = To.copy()
        Ta_n = Ta.copy()
        
        # Update nested grid region with finer resolution
        for i in range(self.offset, self.offset + self.nested_grid_size):
            for j in range(self.offset, self.offset + self.nested_grid_size):
                Q = k * (Ta[i,j] - To[i,j])
                To_n[i,j] += (Q / Co) * dt
                Ta_n[i,j] -= (Q / Ca) * dt
        
        return To_n, Ta_n