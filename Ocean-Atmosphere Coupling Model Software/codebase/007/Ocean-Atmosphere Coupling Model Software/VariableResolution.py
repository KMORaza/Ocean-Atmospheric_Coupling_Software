import numpy as np
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class VariableResolutionGrid:
    def __init__(self, grid_size, coast_factor):
        logging.debug(f"Initializing VariableResolutionGrid with grid_size={grid_size}, coast_factor={coast_factor}")
        try:
            self.grid_size = grid_size
            self.coast_factor = coast_factor
            self.dx = np.ones((grid_size, grid_size))
            self.dy = np.ones((grid_size, grid_size))
            
            # Simulate coastal region with finer resolution
            coast_width = grid_size // 4
            for i in range(grid_size):
                for j in range(grid_size):
                    if i < coast_width or j < coast_width:
                        self.dx[i,j] /= coast_factor
                        self.dy[i,j] /= coast_factor
            logging.debug("VariableResolutionGrid initialization complete")
        except Exception as e:
            logging.error(f"VariableResolutionGrid initialization failed: {str(e)}")
            raise
    
    def get_spatial_steps(self):
        logging.debug("Getting spatial steps")
        try:
            return self.dx, self.dy
        except Exception as e:
            logging.error(f"Get spatial steps failed: {str(e)}")
            raise