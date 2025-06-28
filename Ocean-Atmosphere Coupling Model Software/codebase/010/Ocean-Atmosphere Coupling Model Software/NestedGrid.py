import numpy as np
from scipy.interpolate import RegularGridInterpolator
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class NestedGrid:
    def __init__(self, global_grid_size, nested_grid_size, coast_factor):
        logging.debug(f"Initializing NestedGrid with global_grid_size={global_grid_size}, "
                     f"nested_grid_size={nested_grid_size}, coast_factor={coast_factor}")
        try:
            if nested_grid_size >= global_grid_size:
                raise ValueError("Nested grid size must be smaller than global grid size.")
            if nested_grid_size % 2 != 0:
                raise ValueError("Nested grid size must be even for 2x refinement.")
            
            self.global_grid_size = global_grid_size
            self.nested_grid_size = nested_grid_size
            self.offset = (global_grid_size - nested_grid_size) // 2
            self.nested_dx = 1.0 / (2 * coast_factor)
            
            # Initialize nested grid with 2x resolution
            self.nested_grid_resolution = self.nested_grid_size * 2
            self.nested_ocean = np.zeros((self.nested_grid_resolution, self.nested_grid_resolution))
            self.nested_atm = np.zeros((self.nested_grid_resolution, self.nested_grid_resolution))
            
            logging.debug("NestedGrid initialization complete")
        except Exception as e:
            logging.error(f"NestedGrid initialization failed: {str(e)}")
            raise
    
    def update(self, To, Ta, k, Co, Ca, dt):
        logging.debug("Updating nested grid")
        try:
            # Validate input arrays
            if To.shape != (self.global_grid_size, self.global_grid_size) or \
               Ta.shape != (self.global_grid_size, self.global_grid_size):
                raise ValueError(f"Input arrays must match global grid size ({self.global_grid_size}x{self.global_grid_size}).")
            
            # Clip input arrays to prevent extreme values
            To = np.clip(To, 250, 350)
            Ta = np.clip(Ta, 250, 350)
            
            # Interpolate coarse grid to nested grid
            x_coarse = np.linspace(0, self.global_grid_size-1, self.global_grid_size)
            x_fine = np.linspace(self.offset, self.offset + self.nested_grid_size-1, 
                               self.nested_grid_resolution)
            interp_ocean = RegularGridInterpolator((x_coarse, x_coarse), To, 
                                                method="linear", bounds_error=False, fill_value=0)
            interp_atm = RegularGridInterpolator((x_coarse, x_coarse), Ta, 
                                               method="linear", bounds_error=False, fill_value=0)
            
            X, Y = np.meshgrid(x_fine, x_fine)
            points = np.array([X.flatten(), Y.flatten()]).T
            self.nested_ocean = np.clip(interp_ocean(points).reshape(self.nested_grid_resolution, 
                                                          self.nested_grid_resolution), 250, 350)
            self.nested_atm = np.clip(interp_atm(points).reshape(self.nested_grid_resolution, 
                                                       self.nested_grid_resolution), 250, 350)
            
            # Update nested grid
            Q = k * np.clip(self.nested_atm - self.nested_ocean, -1e5, 1e5)
            self.nested_ocean += np.clip(Q / Co, -1e3, 1e3) * dt
            self.nested_atm -= np.clip(Q / Ca, -1e3, 1e3) * dt
            
            # Clip nested grid temperatures
            self.nested_ocean = np.clip(self.nested_ocean, 250, 350)
            self.nested_atm = np.clip(self.nested_atm, 250, 350)
            
            # Feedback to coarse grid
            To_new = To.copy()
            Ta_new = Ta.copy()
            for i in range(self.nested_grid_size):
                for j in range(self.nested_grid_size):
                    coarse_i = self.offset + i
                    coarse_j = self.offset + j
                    fine_i, fine_j = i * 2, j * 2
                    if (fine_i + 1 < self.nested_grid_resolution and 
                        fine_j + 1 < self.nested_grid_resolution and
                        coarse_i < self.global_grid_size and 
                        coarse_j < self.global_grid_size):
                        To_new[coarse_i, coarse_j] = np.mean(
                            self.nested_ocean[fine_i:fine_i+2, fine_j:fine_j+2]
                        )
                        Ta_new[coarse_i, coarse_j] = np.mean(
                            self.nested_atm[fine_i:fine_i+2, fine_j:fine_j+2]
                        )
            
            logging.debug("Nested grid update complete")
            return To_new, Ta_new
        except Exception as e:
            logging.error(f"Nested grid update failed: {str(e)}")
            raise RuntimeError(f"Nested grid update failed: {str(e)}")