import numpy as np
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class TwoWayCoupling:
    def __init__(self, drag_coeff, wind_speed, precip_rate, evap_rate, solar_forcing, 
                 longwave_coeff, mixing_coeff, co2_transfer_coeff, freshwater_conservation_coeff=1.0, co2_conservation_coeff=1.0):
        logging.debug(f"Initializing TwoWayCoupling with drag_coeff={drag_coeff}, "
                     f"wind_speed={wind_speed}, precip_rate={precip_rate}, "
                     f"evap_rate={evap_rate}, solar_forcing={solar_forcing}, "
                     f"longwave_coeff={longwave_coeff}, mixing_coeff={mixing_coeff}, "
                     f"co2_transfer_coeff={co2_transfer_coeff}, "
                     f"freshwater_conservation_coeff={freshwater_conservation_coeff}, "
                     f"co2_conservation_coeff={co2_conservation_coeff}")
        try:
            self.drag_coeff = drag_coeff
            self.wind_speed = wind_speed
            self.precip_rate = precip_rate
            self.evap_rate = evap_rate
            self.solar_forcing = solar_forcing
            self.longwave_coeff = longwave_coeff
            self.mixing_coeff = mixing_coeff
            self.co2_transfer_coeff = co2_transfer_coeff
            self.freshwater_conservation_coeff = freshwater_conservation_coeff  # New: Scaling factor for mass conservation
            self.co2_conservation_coeff = co2_conservation_coeff  # New: Scaling factor for CO2 conservation
            self.rho_air = 1.225
            self.rho_water = 1025
            self.Cp_air = 1005
            self.Lv = 2.5e6
            self.co2_solubility = 0.03
            logging.debug("TwoWayCoupling initialization complete")
        except Exception as e:
            logging.error(f"TwoWayCoupling initialization failed: {str(e)}")
            raise
    
    def compute_sea_surface_roughness(self, wind_speed, u_ocean, v_ocean):
        logging.debug("Computing sea surface roughness")
        try:
            g = 9.81
            u_star = np.sqrt(self.rho_air * self.drag_coeff * wind_speed**2 / self.rho_water)
            ocean_speed = np.sqrt(u_ocean**2 + v_ocean**2)
            alpha = 0.018
            z0 = alpha * (u_star**2 + 0.1 * ocean_speed**2) / g
            drag_coeff_adj = self.drag_coeff * (1 + 0.1 * np.log10(np.clip(z0, 1e-6, 1e-2)))
            return np.clip(drag_coeff_adj, 1e-4, 1e-2)
        except Exception as e:
            logging.error(f"Sea surface roughness computation failed: {str(e)}")
            raise
    
    def compute_momentum_flux(self, wind_speed, u_ocean, v_ocean):
        logging.debug("Computing momentum flux")
        try:
            drag_coeff = self.compute_sea_surface_roughness(wind_speed, u_ocean, v_ocean)
            tau = self.rho_air * drag_coeff * np.clip(wind_speed**2, 0, 1e3)
            return np.clip(tau, -1e5, 1e5)
        except Exception as e:
            logging.error(f"Momentum flux computation failed: {str(e)}")
            raise
    
    def compute_heat_flux(self, Ta, To, u_ocean, v_ocean):
        logging.debug("Computing heat flux")
        try:
            drag_coeff = self.compute_sea_surface_roughness(self.wind_speed, u_ocean, v_ocean)
            k_sensible = 0.01 * drag_coeff / self.drag_coeff
            Q_sensible = self.rho_air * self.Cp_air * k_sensible * np.clip(Ta - To, -100, 100)
            Q_latent = self.rho_air * self.Lv * self.evap_rate * np.sign(Ta - To)
            Q_total = Q_sensible + Q_latent
            return np.clip(Q_total, -1e6, 1e6)
        except Exception as e:
            logging.error(f"Heat flux computation failed: {str(e)}")
            raise
    
    def compute_freshwater_flux(self, salinity, moisture, ocean_depth=1000.0):
        logging.debug("Computing freshwater flux")
        try:
            # Modified: Ensure mass conservation by balancing precipitation and evaporation
            precip = self.precip_rate * (1 + 0.5 * np.clip(moisture / 0.01, 0, 2))
            evap = self.evap_rate * self.freshwater_conservation_coeff * (1 + 0.1 * np.clip(salinity / 35.0, 0.8, 1.2))
            F = precip - evap  # Net freshwater flux
            # Compute salinity change, ensuring mass conservation
            dS_dt = -salinity * F / (self.rho_water * ocean_depth)  # Normalize by ocean depth
            # Return both salinity change and net flux for conservation tracking
            return np.clip(dS_dt, -1e-3, 1e-3), np.clip(F, -1e-3, 1e-3)
        except Exception as e:
            logging.error(f"Freshwater flux computation failed: {str(e)}")
            raise
    
    def compute_moisture_advection(self, moisture, dx, dy, step, u_atm, v_atm):
        logging.debug("Computing moisture advection")
        try:
            advection = np.zeros_like(moisture)
            for i in range(1, moisture.shape[0]-1):
                for j in range(1, moisture.shape[1]-1):
                    u_ij = u_atm[i,j] if isinstance(u_atm, np.ndarray) else u_atm
                    v_ij = v_atm[i,j] if isinstance(v_atm, np.ndarray) else v_atm
                    adv_x = -u_ij * (moisture[i+1,j] - moisture[i-1,j]) / (2 * dx[i,j])
                    adv_y = -v_ij * (moisture[i,j+1] - moisture[i,j-1]) / (2 * dy[i,j])
                    advection[i,j] = np.clip(adv_x + adv_y, -1e-4, 1e-4)
            return advection
        except Exception as e:
            logging.error(f"Moisture advection computation failed: {str(e)}")
            raise
    
    def compute_turbulent_mixing(self, T, S, dx, dy, wind_speed):
        logging.debug("Computing turbulent mixing")
        try:
            mixing = np.zeros_like(T)
            u_star = np.sqrt(self.rho_air * self.drag_coeff * wind_speed**2 / self.rho_water)
            for i in range(1, T.shape[0]-1):
                for j in range(1, T.shape[1]-1):
                    grad_x = (T[i+1,j] - T[i-1,j]) / (2 * dx[i,j])
                    grad_y = (T[i,j+1] - T[i,j-1]) / (2 * dy[i,j])
                    mixing[i,j] = self.mixing_coeff * u_star * np.clip(grad_x + grad_y, -1e3, 1e3)
            return np.clip(mixing, -1e3, 1e3)
        except Exception as e:
            logging.error(f"Turbulent mixing computation failed: {str(e)}")
            raise
    
    def compute_co2_flux(self, co2_ocean, co2_atm, ocean_depth=1000.0, atm_height=10000.0):
        logging.debug("Computing CO2 flux")
        try:
            # Modified: Ensure CO2 mass conservation between ocean and atmosphere
            pCO2_ocean = co2_ocean / self.co2_solubility
            pCO2_atm = co2_atm
            # Compute flux with conservation scaling
            F_co2 = self.co2_transfer_coeff * self.co2_conservation_coeff * (pCO2_ocean - pCO2_atm)
            # Normalize by layer depths to conserve mass
            F_co2_ocean = F_co2 / ocean_depth  # CO2 flux per unit volume for ocean
            F_co2_atm = -F_co2 / atm_height  # Negative for atmosphere to conserve mass
            return np.clip(F_co2_ocean, -1e-3, 1e-3), np.clip(F_co2_atm, -1e-3, 1e-3)
        except Exception as e:
            logging.error(f"CO2 flux computation failed: {str(e)}")
            raise
    
    def compute_radiative_flux(self, T, co2_atm):
        logging.debug("Computing radiative flux")
        try:
            sigma = 5.67e-8
            T_norm = np.clip(T / 300.0, 0.8, 1.2)
            longwave = self.longwave_coeff * sigma * (T_norm * 300.0)**4
            greenhouse = 0.1 * np.clip(np.log(co2_atm / 400), -1, 1)
            return np.clip(self.solar_forcing - longwave * (1 + greenhouse), -1e6, 1e6)
        except Exception as e:
            logging.error(f"Radiative flux computation failed: {str(e)}")
            raise