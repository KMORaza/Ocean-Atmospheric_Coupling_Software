import numpy as np
import logging
from VariableResolution import VariableResolutionGrid
from NestedGrid import NestedGrid
from AdaptiveMeshRefinement import AdaptiveMeshRefinement

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class OceanAtmosphereModel:
    def __init__(self, ocean_temp, atm_temp, coupling_coeff, heat_capacity_ocean, 
                 heat_capacity_atm, time_step, total_time, grid_size, coast_factor, 
                 use_nested_grid, nested_grid_size, amr_threshold, solar_forcing, 
                 longwave_coeff, adv_velocity):
        logging.debug("Initializing OceanAtmosphereModel")
        try:
            self.To = ocean_temp  # Initial ocean temperature (K)
            self.Ta = atm_temp    # Initial atmosphere temperature (K)
            self.k = coupling_coeff  # Coupling coefficient (W/m^2/K)
            self.Co = heat_capacity_ocean  # Ocean heat capacity (J/m^2/K)
            self.Ca = heat_capacity_atm    # Atmosphere heat capacity (J/m^2/K)
            self.dt = time_step   # Time step (s)
            self.total_time = total_time  # Total simulation time (s)
            self.grid_size = grid_size  # Base grid size (NxN)
            self.solar_forcing = solar_forcing  # Solar forcing (W/m^2)
            self.longwave_coeff = longwave_coeff  # Longwave radiation coefficient
            self.adv_velocity = adv_velocity  # Advection velocity (m/s)
            
            # Initialize grids
            self.grid = VariableResolutionGrid(self.grid_size, coast_factor)
            self.dx, self.dy = self.grid.get_spatial_steps()
            
            self.nested_grid = None
            if use_nested_grid:
                self.nested_grid = NestedGrid(self.grid_size, nested_grid_size, coast_factor)
            
            self.amr = AdaptiveMeshRefinement(self.grid_size, amr_threshold)
            
            # Initialize temperature fields with clipping to prevent extreme values
            self.ocean_temps = np.clip(np.full((self.grid_size, self.grid_size), self.To), 250, 350)
            self.atm_temps = np.clip(np.full((self.grid_size, self.grid_size), self.Ta), 250, 350)
            
            logging.debug("OceanAtmosphereModel initialization complete")
        except Exception as e:
            logging.error(f"OceanAtmosphereModel initialization failed: {str(e)}")
            raise
    
    def step(self, step):
        logging.debug(f"Model step {step}")
        try:
            To = self.ocean_temps.copy()
            Ta = self.atm_temps.copy()
            time = step * self.dt
            
            # Apply AMR
            refinement_mask = self.amr.compute_refinement(To, Ta)
            
            # Update nested grid
            if self.nested_grid:
                To, Ta = self.nested_grid.update(To, Ta, self.k, self.Co, self.Ca, self.dt)
            
            # Compute heat flux: Q = k * (Ta - To)
            Q = self.k * np.clip(Ta - To, -1e5, 1e5)  # Clip to prevent extreme values
            
            # Radiative forcing: solar + longwave (Stefan-Boltzmann approximation)
            sigma = 5.67e-8  # Stefan-Boltzmann constant
            # Normalize temperatures for radiative calculation
            To_norm = np.clip(To / 300.0, 0.8, 1.2)  # Normalize around 300K
            Ta_norm = np.clip(Ta / 300.0, 0.8, 1.2)
            R_ocean = self.solar_forcing - self.longwave_coeff * sigma * (To_norm * 300.0)**4
            R_atm = self.solar_forcing - self.longwave_coeff * sigma * (Ta_norm * 300.0)**4
            
            # Advection
            adv_ocean = self.compute_advection(To, self.dx, self.dy, step)
            adv_atm = self.compute_advection(Ta, self.dx, self.dy, step)
            
            # Diffusion
            diff_ocean = self.compute_diffusion(To, self.dx, self.dy)
            diff_atm = self.compute_diffusion(Ta, self.dx, self.dy)
            
            # Semi-implicit update
            alpha = 0.5  # Crank-Nicolson weighting
            To_new = To + self.dt * (
                alpha * np.clip(Q / self.Co + R_ocean / self.Co - adv_ocean + diff_ocean, -1e3, 1e3) +
                (1 - alpha) * np.clip(Q / self.Co + R_ocean / self.Co, -1e3, 1e3)
            )
            Ta_new = Ta + self.dt * (
                alpha * np.clip(-Q / self.Ca + R_atm / self.Ca - adv_atm + diff_atm, -1e3, 1e3) +
                (1 - alpha) * np.clip(-Q / self.Ca + R_atm / self.Ca, -1e3, 1e3)
            )
            
            # Clip temperatures to realistic range
            self.ocean_temps = np.clip(To_new, 250, 350)
            self.atm_temps = np.clip(Ta_new, 250, 350)
            
            # Apply AMR refinement
            self.ocean_temps = np.where(refinement_mask, self.amr.refine(self.ocean_temps), self.ocean_temps)
            self.atm_temps = np.where(refinement_mask, self.amr.refine(self.atm_temps), self.atm_temps)
            
            logging.debug(f"Model step {step} complete")
            return time, self.ocean_temps, self.atm_temps, refinement_mask
        except Exception as e:
            logging.error(f"Model step {step} failed: {str(e)}")
            raise
    
    def compute_diffusion(self, T, dx, dy):
        logging.debug("Computing diffusion")
        try:
            diffusion = np.zeros_like(T)
            D = 1e-6  # Diffusion coefficient
            for i in range(1, T.shape[0]-1):
                for j in range(1, T.shape[1]-1):
                    diff_x = (T[i+1,j] - 2*T[i,j] + T[i-1,j]) / dx[i,j]**2
                    diff_y = (T[i,j+1] - 2*T[i,j] + T[i,j-1]) / dy[i,j]**2
                    diffusion[i,j] = D * np.clip(diff_x + diff_y, -1e5, 1e5)
            return diffusion
        except Exception as e:
            logging.error(f"Diffusion computation failed: {str(e)}")
            raise
    
    def compute_advection(self, T, dx, dy, step):
        logging.debug("Computing advection")
        try:
            advection = np.zeros_like(T)
            u = self.adv_velocity * np.cos(2 * np.pi * step * self.dt / self.total_time)
            v = self.adv_velocity * np.sin(2 * np.pi * step * self.dt / self.total_time)
            for i in range(1, T.shape[0]-1):
                for j in range(1, T.shape[1]-1):
                    adv_x = -u * (T[i+1,j] - T[i-1,j]) / (2 * dx[i,j])
                    adv_y = -v * (T[i,j+1] - T[i,j-1]) / (2 * dy[i,j])
                    advection[i,j] = np.clip(adv_x + adv_y, -1e5, 1e5)
            return advection
        except Exception as e:
            logging.error(f"Advection computation failed: {str(e)}")
            raise