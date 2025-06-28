import numpy as np

class AdaptiveMeshRefinement:
    def __init__(self, grid_size, threshold):
        self.grid_size = grid_size
        self.threshold = threshold  # Gradient threshold for refinement
        
    def compute_refinement(self, To, Ta):
        # Compute temperature gradients
        grad_To = np.abs(np.gradient(To)[0]) + np.abs(np.gradient(To)[1])
        grad_Ta = np.abs(np.gradient(Ta)[0]) + np.abs(np.gradient(Ta)[1])
        
        # Refine where gradients exceed threshold
        return (grad_To > self.threshold) | (grad_Ta > self.threshold)
    
    def refine(self, T):
        # Simple refinement: average with neighbors for smoother transitions
        T_refined = T.copy()
        for i in range(1, self.grid_size-1):
            for j in range(1, self.grid_size-1):
                T_refined[i,j] = 0.5 * T[i,j] + 0.125 * (
                    T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1]
                )
        return T_refined