import numpy as np
from VariableResolution import VariableResolutionGrid
from NestedGrid import NestedGrid
from AdaptiveMeshRefinement import AdaptiveMeshRefinement

class OceanAtmosphereModel:
    def __init__(self, ocean_temp, atm_temp, coupling_coeff, heat_capacity_ocean, 
                 heat_capacity_atm, time_step, total_time, grid_size, coast_factor, 
                 use_nested_grid, nested_grid_size, amr_threshold):
        self.To = ocean_temp  # Initial ocean temperature (K)
        self.Ta = atm_temp    # Initial atmosphere temperature (K)
        self.k = coupling_coeff  # Coupling coefficient (W/m^2/K)
        self.Co = heat_capacity_ocean  # Ocean heat capacity (J/m^2/K)
        self.Ca = heat_capacity_atm    # Atmosphere heat capacity (J/m^2/K)
        self.dt = time_step   # Time step (s)
        self.total_time = total_time  # Total simulation time (s)
        self.grid_size = grid_size  # Base grid size (NxN)
        
        # Initialize variable resolution grid
        self.grid = VariableResolutionGrid(self.grid_size, coast_factor)
        self.dx, self.dy = self.grid.get_spatial_steps()
        
        # Initialize nested grid if enabled
        self.nested_grid = None
        if use_nested_grid:
            self.nested_grid = NestedGrid(self.grid_size, nested_grid_size, coast_factor)
        
        # Initialize AMR
        self.amr = AdaptiveMeshRefinement(self.grid_size, amr_threshold)
        
        # Initialize temperature fields
        self.ocean_temps = np.full((self.grid_size, self.grid_size), self.To)
        self.atm_temps = np.full((self.grid_size, self.grid_size), self.Ta)
        
    def run(self):
        steps = int(self.total_time / self.dt)
        time = np.linspace(0, self.total_time, steps)
        ocean_temps = np.zeros((steps, self.grid_size, self.grid_size))
        atm_temps = np.zeros((steps, self.grid_size, self.grid_size))
        
        To = self.ocean_temps.copy()
        Ta = self.atm_temps.copy()
        ocean_temps[0] = To
        atm_temps[0] = Ta
        
        for i in range(1, steps):
            # Apply AMR to refine grid where needed
            refinement_mask = self.amr.compute_refinement(To, Ta)
            
            # Update nested grid if enabled
            if self.nested_grid:
                To, Ta = self.nested_grid.update(To, Ta, self.k, self.Co, self.Ca, self.dt)
            
            # Compute heat flux: Q = k * (Ta - To)
            Q = self.k * (Ta - To)
            
            # Temperature changes: dT/dt = Q/C + diffusion
            dTo_dt = Q / self.Co + self.compute_diffusion(To, self.dx, self.dy)
            dTa_dt = -Q / self.Ca + self.compute_diffusion(Ta, self.dx, self.dy)
            
            # Update temperatures using Euler method
            To += dTo_dt * self.dt
            Ta += dTa_dt * self.dt
            
            # Apply refinement mask (simplified: average refined areas)
            To = np.where(refinement_mask, (To + self.amr.refine(To)) / 2, To)
            Ta = np.where(refinement_mask, (Ta + self.amr.refine(Ta)) / 2, Ta)
            
            ocean_temps[i] = To
            atm_temps[i] = Ta
        
        return time, ocean_temps, atm_temps
    
    def compute_diffusion(self, T, dx, dy):
        # Simple 2D diffusion using finite differences
        diffusion = np.zeros_like(T)
        for i in range(1, T.shape[0]-1):
            for j in range(1, T.shape[1]-1):
                diffusion[i,j] = (
                    (T[i+1,j] - 2*T[i,j] + T[i-1,j]) / dx[i,j]**2 +
                    (T[i,j+1] - 2*T[i,j] + T[i,j-1]) / dy[i,j]**2
                ) * 1e-6  # Diffusion coefficient
        return diffusion