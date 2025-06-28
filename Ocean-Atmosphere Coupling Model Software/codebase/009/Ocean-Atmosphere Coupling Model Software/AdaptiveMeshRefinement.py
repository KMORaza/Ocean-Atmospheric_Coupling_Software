import numpy as np
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class AdaptiveMeshRefinement:
    def __init__(self, grid_size, threshold):
        logging.debug(f"Initializing AdaptiveMeshRefinement with grid_size={grid_size}, threshold={threshold}")
        try:
            self.grid_size = grid_size
            self.threshold = threshold
            self.refined_cells = []
            logging.debug("AdaptiveMeshRefinement initialization complete")
        except Exception as e:
            logging.error(f"AdaptiveMeshRefinement initialization failed: {str(e)}")
            raise
    
    def compute_refinement(self, To, Ta):
        logging.debug("Computing refinement")
        try:
            # Clip inputs to prevent extreme values
            To = np.clip(To, 250, 350)
            Ta = np.clip(Ta, 250, 350)
            
            # Compute gradients and vorticity
            grad_To_x, grad_To_y = np.gradient(To)
            grad_Ta_x, grad_Ta_y = np.gradient(Ta)
            grad_mag = np.clip(np.sqrt(np.clip(grad_To_x**2, -1e10, 1e10) + 
                                     np.clip(grad_To_y**2, -1e10, 1e10) + 
                                     np.clip(grad_Ta_x**2, -1e10, 1e10) + 
                                     np.clip(grad_Ta_y**2, -1e10, 1e10)), 0, 1e5)
            
            # Vorticity (simplified for 2D)
            vort_x = np.gradient(grad_To_x)[1]
            vort_y = np.gradient(grad_To_y)[0]
            vort = np.clip(np.abs(vort_x - vort_y), 0, 1e5)
            
            # Identify cells for refinement
            mask = (grad_mag > self.threshold) | (vort > self.threshold / 10)
            self.refined_cells = [(i, j) for i in range(1, self.grid_size-1) 
                                for j in range(1, self.grid_size-1) if mask[i,j]]
            logging.debug("Refinement computation complete")
            return mask
        except Exception as e:
            logging.error(f"Refinement computation failed: {str(e)}")
            raise
    
    def refine(self, T):
        logging.debug("Refining grid")
        try:
            T_refined = T.copy()
            for i, j in self.refined_cells:
                subgrid = np.zeros((2, 2))
                subgrid[0,0] = 0.5 * T[i,j] + 0.125 * (T[i-1,j] + T[i+1,j] + T[i,j-1] + T[i,j+1])
                subgrid[0,1] = 0.5 * T[i,j] + 0.125 * (T[i-1,j] + T[i+1,j] + T[i,j+1] + T[i,j-1])
                subgrid[1,0] = 0.5 * T[i,j] + 0.125 * (T[i+1,j] + T[i-1,j] + T[i,j-1] + T[i,j+1])
                subgrid[1,1] = 0.5 * T[i,j] + 0.125 * (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1])
                T_refined[i,j] = np.clip(np.mean(subgrid), 250, 350)
            logging.debug("Grid refinement complete")
            return T_refined
        except Exception as e:
            logging.error(f"Grid refinement failed: {str(e)}")
            raise