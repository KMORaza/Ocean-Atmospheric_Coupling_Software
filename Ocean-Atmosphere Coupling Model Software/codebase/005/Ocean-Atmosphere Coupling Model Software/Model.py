import numpy as np
import logging
from VariableResolution import VariableResolutionGrid
from NestedGrid import NestedGrid
from AdaptiveMeshRefinement import AdaptiveMeshRefinement
from TwoWayCoupling import TwoWayCoupling

class OceanAtmosphereModel:
    def __init__(self, ocean_temp, atm_temp, coupling_coeff, heat_capacity_ocean, 
                 heat_capacity_atm, time_step, total_time, grid_size, coast_factor, 
                 use_nested_grid, nested_grid_size, amr_threshold, solar_forcing, 
                 longwave_coeff, adv_velocity, drag_coeff, wind_speed, precip_rate, 
                 evap_rate, mixing_coeff, co2_transfer_coeff, freshwater_conservation_coeff=1.0, 
                 co2_conservation_coeff=1.0, ocean_time_scale=1.0, atm_time_scale=0.1):
        
        logging.debug("Initializing OceanAtmosphereModel with time scales: "
                     f"ocean={ocean_time_scale}, atmosphere={atm_time_scale}")
        
        self.To = ocean_temp
        self.Ta = atm_temp
        self.k = coupling_coeff
        self.Co = heat_capacity_ocean
        self.Ca = heat_capacity_atm
        self.dt = time_step
        self.total_time = total_time
        self.grid_size = grid_size
        self.adv_velocity = adv_velocity
        
        # Time scale coupling
        self.ocean_time_scale = ocean_time_scale
        self.atm_time_scale = atm_time_scale
        self.ocean_dt = self.dt * self.ocean_time_scale
        self.atm_dt = self.dt * self.atm_time_scale
        
        self.coupling = TwoWayCoupling(drag_coeff, wind_speed, precip_rate, evap_rate, 
                                     solar_forcing, longwave_coeff, mixing_coeff, 
                                     co2_transfer_coeff, freshwater_conservation_coeff, 
                                     co2_conservation_coeff)
        
        self.grid = VariableResolutionGrid(self.grid_size, coast_factor)
        self.dx, self.dy = self.grid.get_spatial_steps()
        
        self.nested_grid = None
        if use_nested_grid:
            self.nested_grid = NestedGrid(self.grid_size, nested_grid_size, coast_factor)
        
        self.amr = AdaptiveMeshRefinement(self.grid_size, amr_threshold)
        
        # Initialize fields
        self.ocean_temps = np.clip(np.full((self.grid_size, self.grid_size), self.To), 250, 350)
        self.atm_temps = np.clip(np.full((self.grid_size, self.grid_size), self.Ta), 250, 350)
        self.salinity = np.full((self.grid_size, self.grid_size), 35.0)
        self.u_ocean = np.zeros((self.grid_size, self.grid_size))
        self.v_ocean = np.zeros((self.grid_size, self.grid_size))
        self.moisture = np.full((self.grid_size, self.grid_size), 0.01)
        self.co2_ocean = np.full((self.grid_size, self.grid_size), 2.0)
        self.co2_atm = np.full((self.grid_size, self.grid_size), 400.0)
        
        logging.debug("OceanAtmosphereModel initialization complete")

    def step(self, step):
        logging.debug(f"Model step {step}")
        try:
            To = self.ocean_temps.copy()
            Ta = self.atm_temps.copy()
            S = self.salinity.copy()
            u = self.u_ocean.copy()
            v = self.v_ocean.copy()
            M = self.moisture.copy()
            C_o = self.co2_ocean.copy()
            C_a = self.co2_atm.copy()
            time = step * self.dt
            
            refinement_mask = self.amr.compute_refinement(To, Ta)
            
            if self.nested_grid:
                To, Ta = self.nested_grid.update(To, Ta, self.k, self.Co, self.Ca, self.dt)
            
            u_atm = np.full((self.grid_size, self.grid_size), 
                          self.adv_velocity * np.cos(2 * np.pi * step * self.dt / self.total_time))
            v_atm = np.full((self.grid_size, self.grid_size), 
                          self.adv_velocity * np.sin(2 * np.pi * step * self.dt / self.total_time))
            
            Q = self.coupling.compute_heat_flux(Ta, To, u, v)
            R_ocean = self.coupling.compute_radiative_flux(To, C_a)
            R_atm = self.coupling.compute_radiative_flux(Ta, C_a)
            dS_dt, F_freshwater = self.coupling.compute_freshwater_flux(S, M)
            tau = self.coupling.compute_momentum_flux(self.coupling.wind_speed, u, v)
            M_adv = self.coupling.compute_moisture_advection(M, self.dx, self.dy, step, u_atm, v_atm)
            mix_ocean = self.coupling.compute_turbulent_mixing(To, S, self.dx, self.dy, self.coupling.wind_speed)
            mix_atm = self.coupling.compute_turbulent_mixing(Ta, S, self.dx, self.dy, self.coupling.wind_speed)
            F_co2_ocean, F_co2_atm = self.coupling.compute_co2_flux(C_o, C_a)
            
            u_new = u + self.ocean_dt * tau / self.coupling.rho_water
            v_new = v + self.ocean_dt * tau / self.coupling.rho_water
            
            adv_ocean = self.compute_advection(To, self.dx, self.dy, step, u_new, v_new)
            adv_atm = self.compute_advection(Ta, self.dx, self.dy, step, u_atm, v_atm)
            adv_salinity = self.compute_advection(S, self.dx, self.dy, step, u_new, v_new)
            adv_co2_o = self.compute_advection(C_o, self.dx, self.dy, step, u_new, v_new)
            adv_co2_a = self.compute_advection(C_a, self.dx, self.dy, step, u_atm, v_atm)
            
            diff_ocean = self.compute_diffusion(To, self.dx, self.dy)
            diff_atm = self.compute_diffusion(Ta, self.dx, self.dy)
            diff_salinity = self.compute_diffusion(S, self.dx, self.dy)
            
            alpha = 0.5  # Semi-implicit factor
            
            # Apply different time scales
            To_new = To + self.ocean_dt * (
                alpha * np.clip(Q / self.Co + R_ocean / self.Co - adv_ocean + diff_ocean + mix_ocean / self.Co, -1e3, 1e3) +
                (1 - alpha) * np.clip(Q / self.Co + R_ocean / self.Co, -1e3, 1e3)
            )
            Ta_new = Ta + self.atm_dt * (
                alpha * np.clip(-Q / self.Ca + R_atm / self.Ca - adv_atm + diff_atm + mix_atm / self.Ca, -1e3, 1e3) +
                (1 - alpha) * np.clip(-Q / self.Ca + R_atm / self.Ca, -1e3, 1e3)
            )
            
            S_new = S + self.dt * np.clip(dS_dt + diff_salinity + mix_ocean, -1e-2, 1e-2)
            M_new = M + self.dt * np.clip(M_adv - F_freshwater, -1e-4, 1e-4)
            C_o_new = C_o + self.dt * np.clip(F_co2_ocean + adv_co2_o, -1e-2, 1e-2)
            C_a_new = C_a + self.dt * np.clip(F_co2_atm + adv_co2_a, -1e-2, 1e-2)
            
            self.ocean_temps = np.clip(To_new, 250, 350)
            self.atm_temps = np.clip(Ta_new, 250, 350)
            self.salinity = np.clip(S_new, 30, 40)
            self.u_ocean = np.clip(u_new, -10, 10)
            self.v_ocean = np.clip(v_new, -10, 10)
            self.moisture = np.clip(M_new, 0, 0.05)
            self.co2_ocean = np.clip(C_o_new, 0, 10)
            self.co2_atm = np.clip(C_a_new, 200, 1000)
            
            if refinement_mask.any():
                self.ocean_temps = np.where(refinement_mask, self.amr.refine(self.ocean_temps), self.ocean_temps)
                self.atm_temps = np.where(refinement_mask, self.amr.refine(self.atm_temps), self.atm_temps)
                self.salinity = np.where(refinement_mask, self.amr.refine(self.salinity), self.salinity)
            
            logging.debug(f"Model step {step} complete")
            return time, self.ocean_temps, self.atm_temps, refinement_mask
            
        except Exception as e:
            logging.error(f"Model step {step} failed: {str(e)}\n{traceback.format_exc()}")
            raise
    
    def compute_diffusion(self, T, dx, dy):
        logging.debug("Computing diffusion")
        try:
            diffusion = np.zeros_like(T)
            D = 1e-6
            for i in range(1, T.shape[0]-1):
                for j in range(1, T.shape[1]-1):
                    diff_x = (T[i+1,j] - 2*T[i,j] + T[i-1,j]) / dx[i,j]**2
                    diff_y = (T[i,j+1] - 2*T[i,j] + T[i,j-1]) / dy[i,j]**2
                    diffusion[i,j] = D * np.clip(diff_x + diff_y, -1e5, 1e5)
            return diffusion
        except Exception as e:
            logging.error(f"Diffusion computation failed: {str(e)}")
            raise
    
    def compute_advection(self, T, dx, dy, step, u=None, v=None):
        logging.debug("Computing advection")
        try:
            advection = np.zeros_like(T)
            if u is None or v is None:
                u_vel = self.adv_velocity * np.cos(2 * np.pi * step * self.dt / self.total_time)
                v_vel = self.adv_velocity * np.sin(2 * np.pi * step * self.dt / self.total_time)
            else:
                u_vel = u
                v_vel = v
            for i in range(1, T.shape[0]-1):
                for j in range(1, T.shape[1]-1):
                    u_ij = u_vel[i,j] if isinstance(u_vel, np.ndarray) else u_vel
                    v_ij = v_vel[i,j] if isinstance(v_vel, np.ndarray) else v_vel
                    adv_x = -u_ij * (T[i+1,j] - T[i-1,j]) / (2 * dx[i,j])
                    adv_y = -v_ij * (T[i,j+1] - T[i,j-1]) / (2 * dy[i,j])
                    advection[i,j] = np.clip(adv_x + adv_y, -1e5, 1e5)
            return advection
        except Exception as e:
            logging.error(f"Advection computation failed: {str(e)}")
            raise