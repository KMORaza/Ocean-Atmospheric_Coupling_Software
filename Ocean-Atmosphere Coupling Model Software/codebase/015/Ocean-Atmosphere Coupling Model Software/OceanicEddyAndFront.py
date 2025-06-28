import numpy as np
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGroupBox, QFormLayout, QWidget, QMessageBox, QSizePolicy, QCheckBox
from PyQt5.QtGui import QFont, QPainter, QImage, QPen, QBrush, QColor
from PyQt5.QtCore import Qt, QTimer
import logging
from datetime import datetime

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class EddySimulationWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(400, 400)
        self.image = QImage(400, 400, QImage.Format_RGB32)
        self.image.fill(Qt.white)
        self.vorticity = None
        self.grid_size = 100
        self.time = 0.0
        self.eddy_centers = []
        self.eddy_strengths = []
        self.eddy_radii = []
        self.fine_grids = []
        self.fine_grid_size = 20
        self.refinement_factor = 4
        self.vertical_velocity = None
        self.pressure_perturbation = None
        self.density = None
        self.velocity_shear = None
        self.buoyancy = None
        self.potential_vorticity = None
        self.enable_nonhydrostatic = False
        self.enable_baroclinic = False
        self.enable_barotropic = False
        self.enable_mli = False
        self.enable_pv_frontogenesis = False
        self.grid_spacing_km = 1.0

    def initialize_field(self, grid_size, eddy_strength, eddy_radius, num_eddies, refinement_threshold, 
                        refinement_factor, coriolis_param, rossby_number, enable_nonhydrostatic, 
                        grid_spacing_km, enable_baroclinic, enable_barotropic, baroclinic_shear, 
                        stratification_param, enable_mli, mixed_layer_depth, buoyancy_gradient, 
                        enable_pv_frontogenesis, pv_gradient_factor):
        self.grid_size = grid_size
        self.refinement_factor = refinement_factor
        self.enable_nonhydrostatic = enable_nonhydrostatic
        self.grid_spacing_km = grid_spacing_km
        self.enable_baroclinic = enable_baroclinic
        self.enable_barotropic = enable_barotropic
        self.enable_mli = enable_mli
        self.enable_pv_frontogenesis = enable_pv_frontogenesis
        x = np.linspace(-1, 1, grid_size)
        y = np.linspace(-1, 1, grid_size)
        X, Y = np.meshgrid(x, y)
        self.vorticity = np.zeros((grid_size, grid_size))
        
        # Initialize non-hydrostatic fields
        if self.enable_nonhydrostatic and grid_spacing_km < 1.0:
            self.vertical_velocity = np.zeros((grid_size, grid_size))
            self.pressure_perturbation = np.zeros((grid_size, grid_size))
        
        # Initialize density field for baroclinic instability
        if self.enable_baroclinic:
            self.density = np.ones((grid_size, grid_size)) * 1025.0
            self.density += baroclinic_shear * Y * stratification_param
            self.density = np.clip(self.density, 1020.0, 1030.0)
        
        # Initialize velocity shear for barotropic instability
        if self.enable_barotropic:
            self.velocity_shear = np.zeros((grid_size, grid_size))
            self.velocity_shear += baroclinic_shear * Y
            self.velocity_shear = np.clip(self.velocity_shear, -0.1, 0.1)
        
        # Initialize buoyancy field for mixed-layer instability
        if self.enable_mli and grid_spacing_km < 10.0:
            self.buoyancy = np.zeros((grid_size, grid_size))
            self.buoyancy += buoyancy_gradient * Y
            self.buoyancy += 0.01 * np.sin(2 * np.pi * X / 0.2)
            self.buoyancy = np.clip(self.buoyancy, -0.1, 0.1)
        
        # Initialize potential vorticity field
        if self.enable_pv_frontogenesis:
            self.potential_vorticity = np.zeros((grid_size, grid_size))
            if self.enable_mli and grid_spacing_km < 10.0:
                db_dz = self.buoyancy / mixed_layer_depth
                self.potential_vorticity = self.vorticity + coriolis_param * db_dz / stratification_param
            else:
                self.potential_vorticity = self.vorticity
            self.potential_vorticity = np.clip(self.potential_vorticity, -1.0, 1.0)
        
        # Initialize eddies
        self.eddy_centers = []
        self.eddy_strengths = []
        self.eddy_radii = []
        self.fine_grids = []
        np.random.seed(42)
        for _ in range(num_eddies):
            x0 = np.random.uniform(-0.8, 0.8)
            y0 = np.random.uniform(-0.8, 0.8)
            strength = eddy_strength * np.random.uniform(0.8, 1.2)
            radius = eddy_radius * np.random.uniform(0.8, 1.2)
            self.eddy_centers.append([x0, y0])
            self.eddy_strengths.append(strength)
            self.eddy_radii.append(radius)
            r = np.sqrt((X - x0)**2 + (Y - y0)**2)
            vortex = np.where(r < radius, strength * (r / radius), strength * (radius / r))
            ageostrophic = strength * np.exp(-r**2 / (2 * radius**2)) * (1 - rossby_number * np.tanh(r / radius))
            self.vorticity += vortex * np.exp(-r**2 / (2 * radius**2))
            
            if self.enable_nonhydrostatic and grid_spacing_km < 1.0:
                w = strength * 0.1 * np.exp(-r**2 / (radius**2)) * np.sin(np.arctan2(Y - y0, X - x0))
                p = strength * 0.05 * np.exp(-r**2 / (radius**2))
                self.vertical_velocity += w
                self.pressure_perturbation += p
            
            if self.enable_baroclinic:
                self.density += strength * 0.01 * np.exp(-r**2 / (2 * radius**2))
                self.density = np.clip(self.density, 1020.0, 1030.0)
            
            if self.enable_barotropic:
                self.velocity_shear += strength * 0.05 * np.exp(-r**2 / (2 * radius**2))
                self.velocity_shear = np.clip(self.velocity_shear, -0.1, 0.1)
            
            if self.enable_mli and grid_spacing_km < 10.0:
                self.buoyancy += strength * 0.02 * np.exp(-r**2 / (2 * radius**2))
                self.buoyancy = np.clip(self.buoyancy, -0.1, 0.1)
            
            if self.enable_pv_frontogenesis:
                pv_perturbation = strength * pv_gradient_factor * np.exp(-r**2 / (2 * radius**2))
                self.potential_vorticity += pv_perturbation
                self.potential_vorticity = np.clip(self.potential_vorticity, -1.0, 1.0)
        
        # Initialize fine grid patches
        self.fine_grid_size = max(10, grid_size // 5)
        for (x0, y0), strength, radius in zip(self.eddy_centers, self.eddy_strengths, self.eddy_radii):
            if abs(strength) > refinement_threshold:
                fine_grid = np.zeros((self.fine_grid_size * refinement_factor, self.fine_grid_size * refinement_factor))
                fine_x = np.linspace(x0 - radius, x0 + radius, self.fine_grid_size * refinement_factor)
                fine_y = np.linspace(y0 - radius, y0 + radius, self.fine_grid_size * refinement_factor)
                fine_X, fine_Y = np.meshgrid(fine_x, fine_y)
                r = np.sqrt((fine_X - x0)**2 + (fine_Y - y0)**2)
                vortex = np.where(r < radius, strength * (r / radius), strength * (radius / r))
                fine_vorticity = vortex * np.exp(-r**2 / (2 * radius**2))
                fine_grid_dict = {'vorticity': fine_vorticity, 'center': [x0, y0], 'radius': radius}
                
                if self.enable_nonhydrostatic and grid_spacing_km < 1.0:
                    fine_w = strength * 0.1 * np.exp(-r**2 / (radius**2)) * np.sin(np.arctan2(fine_Y - y0, fine_X - x0))
                    fine_p = strength * 0.05 * np.exp(-r**2 / (radius**2))
                    fine_grid_dict['vertical_velocity'] = fine_w
                    fine_grid_dict['pressure_perturbation'] = fine_p
                
                if self.enable_baroclinic:
                    fine_density = np.ones_like(fine_vorticity) * 1025.0 + strength * 0.01 * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['density'] = np.clip(fine_density, 1020.0, 1030.0)
                
                if self.enable_barotropic:
                    fine_shear = strength * 0.05 * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['velocity_shear'] = np.clip(fine_shear, -0.1, 0.1)
                
                if self.enable_mli and grid_spacing_km < 10.0:
                    fine_buoyancy = strength * 0.02 * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['buoyancy'] = np.clip(fine_buoyancy, -0.1, 0.1)
                
                if self.enable_pv_frontogenesis:
                    fine_pv = fine_vorticity
                    if self.enable_mli and grid_spacing_km < 10.0:
                        db_dz = fine_buoyancy / mixed_layer_depth
                        fine_pv += coriolis_param * db_dz / stratification_param
                    fine_pv += strength * pv_gradient_factor * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['potential_vorticity'] = np.clip(fine_pv, -1.0, 1.0)
                
                self.fine_grids.append(fine_grid_dict)
        
        self.update_visualization()

    def compute_induced_velocity(self, x, y, eddy_centers, eddy_strengths, coriolis_param, rossby_number):
        """Compute velocity at (x, y) with ageostrophic effects."""
        u, v = 0.0, 0.0
        for (x0, y0), strength in zip(eddy_centers, eddy_strengths):
            r_squared = (x - x0)**2 + (y - y0)**2
            if r_squared > 1e-6:
                u_geo = -strength * (y - y0) / (2 * np.pi * r_squared)
                v_geo = strength * (x - x0) / (2 * np.pi * r_squared)
                u_ageo = -rossby_number * v_geo
                v_ageo = rossby_number * u_geo
                u += u_geo + coriolis_param * u_ageo
                v += v_geo + coriolis_param * v_ageo
        return u, v

    def compute_nonhydrostatic_effects(self, vorticity, grid_spacing_km):
        """Compute non-hydrostatic vertical velocity and pressure perturbation."""
        if not self.enable_nonhydrostatic or grid_spacing_km >= 1.0:
            return np.zeros_like(vorticity), np.zeros_like(vorticity)
        
        w = np.zeros_like(vorticity)
        p = np.zeros_like(vorticity)
        dx = 2.0 / self.grid_size
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                u_x = (np.roll(vorticity, -1, axis=1)[i, j] - np.roll(vorticity, 1, axis=1)[i, j]) / (2 * dx)
                v_y = (np.roll(vorticity, -1, axis=0)[i, j] - np.roll(vorticity, 1, axis=0)[i, j]) / (2 * dx)
                w[i, j] = -0.1 * (u_x + v_y)
                p[i, j] = -0.05 * vorticity[i, j]
        return w, p

    def compute_baroclinic_effects(self, density, dx):
        """Compute baroclinic vorticity tendency."""
        if not self.enable_baroclinic:
            return np.zeros_like(density)
        
        grad_rho_x = (np.roll(density, -1, axis=1) - np.roll(density, 1, axis=1)) / (2 * dx)
        grad_rho_y = (np.roll(density, -1, axis=0) - np.roll(density, 1, axis=0)) / (2 * dx)
        vorticity_tendency = -0.1 * (grad_rho_x + grad_rho_y)
        vorticity_tendency = np.clip(vorticity_tendency, -0.1, 0.1)
        return vorticity_tendency

    def compute_barotropic_effects(self, velocity_shear, dx):
        """Compute barotropic vorticity tendency."""
        if not self.enable_barotropic:
            return np.zeros_like(velocity_shear)
        
        grad_shear_x = (np.roll(velocity_shear, -1, axis=1) - np.roll(velocity_shear, 1, axis=1)) / (2 * dx)
        grad_shear_y = (np.roll(velocity_shear, -1, axis=0) - np.roll(velocity_shear, 1, axis=0)) / (2 * dx)
        vorticity_tendency = 0.05 * (grad_shear_x - grad_shear_y)
        vorticity_tendency = np.clip(vorticity_tendency, -0.1, 0.1)
        return vorticity_tendency

    def compute_mli_effects(self, buoyancy, dx, mixed_layer_depth):
        """Compute mixed-layer instability vorticity tendency."""
        if not self.enable_mli or self.grid_spacing_km >= 10.0:
            return np.zeros_like(buoyancy)
        
        grad_b_x = (np.roll(buoyancy, -1, axis=1) - np.roll(buoyancy, 1, axis=1)) / (2 * dx)
        grad_b_y = (np.roll(buoyancy, -1, axis=0) - np.roll(buoyancy, 1, axis=0)) / (2 * dx)
        vorticity_tendency = 0.1 * mixed_layer_depth * (grad_b_x - grad_b_y)
        vorticity_tendency = np.clip(vorticity_tendency, -0.1, 0.1)
        return vorticity_tendency

    def compute_pv_frontogenesis_effects(self, potential_vorticity, dx, strain_rate):
        """Compute potential vorticity frontogenesis tendency."""
        if not self.enable_pv_frontogenesis:
            return np.zeros_like(potential_vorticity)
        
        grad_pv_x = (np.roll(potential_vorticity, -1, axis=1) - np.roll(potential_vorticity, 1, axis=1)) / (2 * dx)
        grad_pv_y = (np.roll(potential_vorticity, -1, axis=0) - np.roll(potential_vorticity, 1, axis=0)) / (2 * dx)
        frontogenesis_tendency = -strain_rate * (grad_pv_x**2 + grad_pv_y**2)
        frontogenesis_tendency = np.clip(frontogenesis_tendency, -0.1, 0.1)
        return frontogenesis_tendency

    def get_simulation_data(self):
        """Collect simulation data for export."""
        data = []
        data.append(f"Simulation Time: {self.time:.2f} s")
        data.append(f"Grid Size: {self.grid_size}")
        data.append(f"Number of Eddies: {len(self.eddy_centers)}")
        data.append(f"Refinement Factor: {self.refinement_factor}")
        data.append(f"Non-Hydrostatic: {self.enable_nonhydrostatic}")
        data.append(f"Baroclinic Instability: {self.enable_baroclinic}")
        data.append(f"Barotropic Instability: {self.enable_barotropic}")
        data.append(f"Mixed-Layer Instability: {self.enable_mli}")
        data.append(f"PV Frontogenesis: {self.enable_pv_frontogenesis}")
        data.append(f"Grid Spacing: {self.grid_spacing_km:.2f} km")
        data.append("\nEddy Information:")
        for i, (center, strength, radius) in enumerate(zip(self.eddy_centers, self.eddy_strengths, self.eddy_radii)):
            data.append(f"Eddy {i+1}: Center=({center[0]:.3f}, {center[1]:.3f}), Strength={strength:.3f}, Radius={radius:.3f}")
        
        data.append("\nCoarse Grid Vorticity Statistics:")
        if self.vorticity is not None:
            data.append(f"Min Vorticity: {np.min(self.vorticity):.3f}")
            data.append(f"Max Vorticity: {np.max(self.vorticity):.3f}")
            data.append(f"Mean Vorticity: {np.mean(self.vorticity):.3f}")
            data.append(f"Std Vorticity: {np.std(self.vorticity):.3f}")
            grad_x = (np.roll(self.vorticity, -1, axis=1) - np.roll(self.vorticity, 1, axis=1)) / (2 * 2/self.grid_size)
            grad_y = (np.roll(self.vorticity, -1, axis=0) - np.roll(self.vorticity, 1, axis=0)) / (2 * 2/self.grid_size)
            grad_magnitude = np.sqrt(grad_x**2 + grad_y**2 + 1e-10)
            data.append(f"Max Vorticity Gradient: {np.max(grad_magnitude):.3f}")
        
        if self.enable_nonhydrostatic and self.grid_spacing_km < 1.0:
            data.append("\nCoarse Grid Non-Hydrostatic Statistics:")
            data.append(f"Min Vertical Velocity: {np.min(self.vertical_velocity):.3f}")
            data.append(f"Max Vertical Velocity: {np.max(self.vertical_velocity):.3f}")
            data.append(f"Mean Vertical Velocity: {np.mean(self.vertical_velocity):.3f}")
            data.append(f"Min Pressure Perturbation: {np.min(self.pressure_perturbation):.3f}")
            data.append(f"Max Pressure Perturbation: {np.max(self.pressure_perturbation):.3f}")
        
        if self.enable_baroclinic:
            data.append("\nCoarse Grid Density Statistics:")
            data.append(f"Min Density: {np.min(self.density):.3f}")
            data.append(f"Max Density: {np.max(self.density):.3f}")
            data.append(f"Mean Density: {np.mean(self.density):.3f}")
        
        if self.enable_barotropic:
            data.append("\nCoarse Grid Velocity Shear Statistics:")
            data.append(f"Min Shear: {np.min(self.velocity_shear):.3f}")
            data.append(f"Max Shear: {np.max(self.velocity_shear):.3f}")
            data.append(f"Mean Shear: {np.mean(self.velocity_shear):.3f}")
        
        if self.enable_mli and self.grid_spacing_km < 10.0:
            data.append("\nCoarse Grid Buoyancy Statistics:")
            data.append(f"Min Buoyancy: {np.min(self.buoyancy):.3f}")
            data.append(f"Max Buoyancy: {np.max(self.buoyancy):.3f}")
            data.append(f"Mean Buoyancy: {np.mean(self.buoyancy):.3f}")
        
        if self.enable_pv_frontogenesis:
            data.append("\nCoarse Grid Potential Vorticity Statistics:")
            data.append(f"Min PV: {np.min(self.potential_vorticity):.3f}")
            data.append(f"Max PV: {np.max(self.potential_vorticity):.3f}")
            data.append(f"Mean PV: {np.mean(self.potential_vorticity):.3f}")
        
        data.append("\nFine Grid Statistics:")
        for i, fine_grid_info in enumerate(self.fine_grids):
            fine_vorticity = fine_grid_info['vorticity']
            x0, y0 = fine_grid_info['center']
            radius = fine_grid_info['radius']
            data.append(f"Fine Grid {i+1}: Center=({x0:.3f}, {y0:.3f}), Radius={radius:.3f}")
            data.append(f"  Min Vorticity: {np.min(fine_vorticity):.3f}")
            data.append(f"  Max Vorticity: {np.max(fine_vorticity):.3f}")
            data.append(f"  Mean Vorticity: {np.mean(fine_vorticity):.3f}")
            data.append(f"  Std Vorticity: {np.std(fine_vorticity):.3f}")
            grad_fine_x = (np.roll(fine_vorticity, -1, axis=1) - np.roll(fine_vorticity, 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
            grad_fine_y = (np.roll(fine_vorticity, -1, axis=0) - np.roll(fine_vorticity, 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
            grad_fine = np.sqrt(grad_fine_x**2 + grad_fine_y**2 + 1e-10)
            data.append(f"  Max Vorticity Gradient: {np.max(grad_fine):.3f}")
            if self.enable_nonhydrostatic and self.grid_spacing_km < 1.0:
                fine_w = fine_grid_info['vertical_velocity']
                fine_p = fine_grid_info['pressure_perturbation']
                data.append(f"  Min Vertical Velocity: {np.min(fine_w):.3f}")
                data.append(f"  Max Vertical Velocity: {np.max(fine_w):.3f}")
                data.append(f"  Min Pressure Perturbation: {np.min(fine_p):.3f}")
                data.append(f"  Max Pressure Perturbation: {np.max(fine_p):.3f}")
            if self.enable_baroclinic:
                fine_density = fine_grid_info['density']
                data.append(f"  Min Density: {np.min(fine_density):.3f}")
                data.append(f"  Max Density: {np.max(fine_density):.3f}")
            if self.enable_barotropic:
                fine_shear = fine_grid_info['velocity_shear']
                data.append(f"  Min Shear: {np.min(fine_shear):.3f}")
                data.append(f"  Max Shear: {np.max(fine_shear):.3f}")
            if self.enable_mli and self.grid_spacing_km < 10.0:
                fine_buoyancy = fine_grid_info['buoyancy']
                data.append(f"  Min Buoyancy: {np.min(fine_buoyancy):.3f}")
                data.append(f"  Max Buoyancy: {np.max(fine_buoyancy):.3f}")
            if self.enable_pv_frontogenesis:
                fine_pv = fine_grid_info['potential_vorticity']
                data.append(f"  Min PV: {np.min(fine_pv):.3f}")
                data.append(f"  Max PV: {np.max(fine_pv):.3f}")
        
        return data

    def update_simulation(self, vorticity_diffusion, rotation_rate, background_flow, decay_rate, 
                         refinement_threshold, strain_rate, coriolis_param, rossby_number, 
                         enable_nonhydrostatic, grid_spacing_km, enable_baroclinic, enable_barotropic, 
                         baroclinic_shear, stratification_param, enable_mli, mixed_layer_depth, 
                         buoyancy_gradient, enable_pv_frontogenesis, pv_gradient_factor):
        if self.vorticity is None:
            return
        dt = 0.1
        max_velocity = max(abs(background_flow), 0.1)
        dx = 2.0 / self.grid_size
        if dt * max_velocity / dx > 0.5:
            dt = 0.5 * dx / max_velocity
        self.time += dt
        x = np.linspace(-1, 1, self.grid_size)
        X, Y = np.meshgrid(x, x)
        dx = 2.0 / self.grid_size
        
        # Update eddy positions
        new_centers = []
        for i, (x0, y0) in enumerate(self.eddy_centers):
            u_ind, v_ind = self.compute_induced_velocity(x0, y0, 
                [c for j, c in enumerate(self.eddy_centers) if j != i],
                [s for j, s in enumerate(self.eddy_strengths) if j != i],
                coriolis_param, rossby_number)
            x0 += (background_flow + u_ind) * dt
            y0 += v_ind * dt
            if x0 > 1:
                x0 -= 2
            new_centers.append([x0, y0])
        self.eddy_centers = new_centers
        
        # Update coarse grid vorticity
        self.vorticity *= np.exp(-decay_rate * dt)
        new_vorticity = np.zeros_like(self.vorticity)
        for (x0, y0), strength, radius in zip(self.eddy_centers, self.eddy_strengths, self.eddy_radii):
            r = np.sqrt((X - x0)**2 + (Y - y0)**2)
            vortex = np.where(r < radius, strength * (r / radius), strength * (radius / r))
            new_vorticity += vortex * np.exp(-r**2 / (2 * radius**2)) * (1 - rossby_number * np.tanh(r / radius))
        
        grad_x = (np.roll(new_vorticity, -1, axis=1) - np.roll(new_vorticity, 1, axis=1)) / (2 * dx)
        grad_y = (np.roll(new_vorticity, -1, axis=0) - np.roll(new_vorticity, 1, axis=0)) / (2 * dx)
        self.vorticity = new_vorticity + vorticity_diffusion * dt * (
            np.roll(new_vorticity, 1, axis=0) + 
            np.roll(new_vorticity, -1, axis=0) +
            np.roll(new_vorticity, 1, axis=1) +
            np.roll(new_vorticity, -1, axis=1) - 
            4 * new_vorticity
        ) - strain_rate * dt * (X * grad_x + Y * grad_y)
        self.vorticity += background_flow * (Y * 0.1)
        self.vorticity = np.clip(self.vorticity, -10.0, 10.0)
        
        # Apply instabilities
        if self.enable_baroclinic:
            vorticity_tendency = self.compute_baroclinic_effects(self.density, dx)
            self.vorticity += dt * vorticity_tendency
            grad_rho_x = (np.roll(self.density, -1, axis=1) - np.roll(self.density, 1, axis=1)) / (2 * dx)
            grad_rho_y = (np.roll(self.density, -1, axis=0) - np.roll(self.density, 1, axis=0)) / (2 * dx)
            self.density -= dt * (grad_x * grad_rho_x + grad_y * grad_rho_y)
            self.density = np.clip(self.density, 1020.0, 1030.0)
        
        if self.enable_barotropic:
            vorticity_tendency = self.compute_barotropic_effects(self.velocity_shear, dx)
            self.vorticity += dt * vorticity_tendency
            grad_shear_x = (np.roll(self.velocity_shear, -1, axis=1) - np.roll(self.velocity_shear, 1, axis=1)) / (2 * dx)
            grad_shear_y = (np.roll(self.velocity_shear, -1, axis=0) - np.roll(self.velocity_shear, 1, axis=0)) / (2 * dx)
            self.velocity_shear -= dt * (grad_x * grad_shear_x + grad_y * grad_shear_y)
            self.velocity_shear = np.clip(self.velocity_shear, -0.1, 0.1)
        
        if self.enable_mli and self.grid_spacing_km < 10.0:
            vorticity_tendency = self.compute_mli_effects(self.buoyancy, dx, mixed_layer_depth)
            self.vorticity += dt * vorticity_tendency
            grad_b_x = (np.roll(self.buoyancy, -1, axis=1) - np.roll(self.buoyancy, 1, axis=1)) / (2 * dx)
            grad_b_y = (np.roll(self.buoyancy, -1, axis=0) - np.roll(self.buoyancy, 1, axis=0)) / (2 * dx)
            self.buoyancy -= dt * (grad_x * grad_b_x + grad_y * grad_b_y)
            self.buoyancy = np.clip(self.buoyancy, -0.1, 0.1)
        
        if self.enable_pv_frontogenesis:
            pv_tendency = self.compute_pv_frontogenesis_effects(self.potential_vorticity, dx, strain_rate)
            self.potential_vorticity += dt * pv_tendency
            self.potential_vorticity = np.clip(self.potential_vorticity, -1.0, 1.0)
            vorticity_tendency = 0.05 * pv_tendency
            self.vorticity += dt * vorticity_tendency
            grad_pv_x = (np.roll(self.potential_vorticity, -1, axis=1) - np.roll(self.potential_vorticity, 1, axis=1)) / (2 * dx)
            grad_pv_y = (np.roll(self.potential_vorticity, -1, axis=0) - np.roll(self.potential_vorticity, 1, axis=0)) / (2 * dx)
            self.potential_vorticity -= dt * (grad_x * grad_pv_x + grad_y * grad_pv_y)
            self.potential_vorticity = np.clip(self.potential_vorticity, -1.0, 1.0)
        
        if self.enable_nonhydrostatic and grid_spacing_km < 1.0:
            self.vertical_velocity, self.pressure_perturbation = self.compute_nonhydrostatic_effects(self.vorticity, grid_spacing_km)
        
        # Update fine grid patches
        self.fine_grids = []
        for (x0, y0), strength, radius in zip(self.eddy_centers, self.eddy_strengths, self.eddy_radii):
            if abs(strength) > refinement_threshold:
                fine_grid = np.zeros((self.fine_grid_size * self.refinement_factor, self.fine_grid_size * self.refinement_factor))
                fine_x = np.linspace(x0 - radius, x0 + radius, self.fine_grid_size * self.refinement_factor)
                fine_y = np.linspace(y0 - radius, y0 + radius, self.fine_grid_size * self.refinement_factor)
                fine_X, fine_Y = np.meshgrid(fine_x, fine_y)
                r = np.sqrt((fine_X - x0)**2 + (fine_Y - y0)**2)
                vortex = np.where(r < radius, strength * (r / radius), strength * (radius / r))
                fine_vorticity = vortex * np.exp(-r**2 / (2 * radius**2)) * (1 - rossby_number * np.tanh(r / radius))
                fine_grid_dict = {'vorticity': fine_vorticity, 'center': [x0, y0], 'radius': radius}
                
                if self.enable_nonhydrostatic and grid_spacing_km < 1.0:
                    fine_w, fine_p = self.compute_nonhydrostatic_effects(fine_vorticity, grid_spacing_km / self.refinement_factor)
                    fine_grid_dict['vertical_velocity'] = fine_w
                    fine_grid_dict['pressure_perturbation'] = fine_p
                
                if self.enable_baroclinic:
                    fine_density = np.ones_like(fine_vorticity) * 1025.0 + strength * 0.01 * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['density'] = np.clip(fine_density, 1020.0, 1030.0)
                
                if self.enable_barotropic:
                    fine_shear = strength * 0.05 * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['velocity_shear'] = np.clip(fine_shear, -0.1, 0.1)
                
                if self.enable_mli and grid_spacing_km < 10.0:
                    fine_buoyancy = strength * 0.02 * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['buoyancy'] = np.clip(fine_buoyancy, -0.1, 0.1)
                
                if self.enable_pv_frontogenesis:
                    fine_pv = fine_vorticity
                    if self.enable_mli and grid_spacing_km < 10.0:
                        db_dz = fine_buoyancy / mixed_layer_depth
                        fine_pv += coriolis_param * db_dz / stratification_param
                    fine_pv += strength * pv_gradient_factor * np.exp(-r**2 / (2 * radius**2))
                    fine_grid_dict['potential_vorticity'] = np.clip(fine_pv, -1.0, 1.0)
                
                grad_fine_x = (np.roll(fine_vorticity, -1, axis=1) - np.roll(fine_vorticity, 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                grad_fine_y = (np.roll(fine_vorticity, -1, axis=0) - np.roll(fine_vorticity, 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                fine_grid = fine_vorticity + (vorticity_diffusion + strain_rate * 0.5) * dt * (
                    np.roll(fine_vorticity, 1, axis=0) + 
                    np.roll(fine_vorticity, -1, axis=0) +
                    np.roll(fine_vorticity, 1, axis=1) +
                    np.roll(fine_vorticity, -1, axis=1) - 
                    4 * fine_vorticity
                ) - strain_rate * dt * (fine_X * grad_fine_x + fine_Y * grad_fine_y)
                fine_grid = np.clip(fine_grid, -10.0, 10.0)
                
                if self.enable_baroclinic:
                    fine_vorticity_tendency = self.compute_baroclinic_effects(fine_grid_dict['density'], 2.0 / (self.fine_grid_size * self.refinement_factor))
                    fine_grid += dt * fine_vorticity_tendency
                    grad_rho_x = (np.roll(fine_grid_dict['density'], -1, axis=1) - np.roll(fine_grid_dict['density'], 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    grad_rho_y = (np.roll(fine_grid_dict['density'], -1, axis=0) - np.roll(fine_grid_dict['density'], 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    fine_grid_dict['density'] -= dt * (grad_fine_x * grad_rho_x + grad_fine_y * grad_rho_y)
                    fine_grid_dict['density'] = np.clip(fine_grid_dict['density'], 1020.0, 1030.0)
                
                if self.enable_barotropic:
                    fine_vorticity_tendency = self.compute_barotropic_effects(fine_grid_dict['velocity_shear'], 2.0 / (self.fine_grid_size * self.refinement_factor))
                    fine_grid += dt * fine_vorticity_tendency
                    grad_shear_x = (np.roll(fine_grid_dict['velocity_shear'], -1, axis=1) - np.roll(fine_grid_dict['velocity_shear'], 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    grad_shear_y = (np.roll(fine_grid_dict['velocity_shear'], -1, axis=0) - np.roll(fine_grid_dict['velocity_shear'], 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    fine_grid_dict['velocity_shear'] -= dt * (grad_fine_x * grad_shear_x + grad_fine_y * grad_shear_y)
                    fine_grid_dict['velocity_shear'] = np.clip(fine_grid_dict['velocity_shear'], -0.1, 0.1)
                
                if self.enable_mli and grid_spacing_km < 10.0:
                    fine_vorticity_tendency = self.compute_mli_effects(fine_grid_dict['buoyancy'], 2.0 / (self.fine_grid_size * self.refinement_factor), mixed_layer_depth)
                    fine_grid += dt * fine_vorticity_tendency
                    grad_b_x = (np.roll(fine_grid_dict['buoyancy'], -1, axis=1) - np.roll(fine_grid_dict['buoyancy'], 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    grad_b_y = (np.roll(fine_grid_dict['buoyancy'], -1, axis=0) - np.roll(fine_grid_dict['buoyancy'], 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    fine_grid_dict['buoyancy'] -= dt * (grad_fine_x * grad_b_x + grad_fine_y * grad_b_y)
                    fine_grid_dict['buoyancy'] = np.clip(fine_grid_dict['buoyancy'], -0.1, 0.1)
                
                if self.enable_pv_frontogenesis:
                    fine_pv_tendency = self.compute_pv_frontogenesis_effects(fine_grid_dict['potential_vorticity'], 2.0 / (self.fine_grid_size * self.refinement_factor), strain_rate)
                    fine_grid_dict['potential_vorticity'] += dt * fine_pv_tendency
                    fine_grid_dict['potential_vorticity'] = np.clip(fine_grid_dict['potential_vorticity'], -1.0, 1.0)
                    fine_vorticity_tendency = 0.05 * fine_pv_tendency
                    fine_grid += dt * fine_vorticity_tendency
                    grad_pv_x = (np.roll(fine_grid_dict['potential_vorticity'], -1, axis=1) - np.roll(fine_grid_dict['potential_vorticity'], 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    grad_pv_y = (np.roll(fine_grid_dict['potential_vorticity'], -1, axis=0) - np.roll(fine_grid_dict['potential_vorticity'], 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                    fine_grid_dict['potential_vorticity'] -= dt * (grad_fine_x * grad_pv_x + grad_fine_y * grad_pv_y)
                    fine_grid_dict['potential_vorticity'] = np.clip(fine_grid_dict['potential_vorticity'], -1.0, 1.0)
                
                fine_grid_dict['vorticity'] = fine_grid
                self.fine_grids.append(fine_grid_dict)
        
        # Interpolate fine grids to coarse grid
        for fine_grid_info in self.fine_grids:
            fine_vorticity = fine_grid_info['vorticity']
            x0, y0 = fine_grid_info['center']
            radius = fine_grid_info['radius']
            i_start = int((x0 - radius + 1) * self.grid_size / 2)
            i_end = i_start + self.fine_grid_size
            j_start = int((y0 - radius + 1) * self.grid_size / 2)
            j_end = j_start + self.fine_grid_size
            if i_start >= 0 and i_end <= self.grid_size and j_start >= 0 and j_end <= self.grid_size:
                coarse_patch = fine_vorticity[::self.refinement_factor, ::self.refinement_factor]
                patch_shape = coarse_patch.shape
                self.vorticity[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_patch
                if self.enable_nonhydrostatic and grid_spacing_km < 1.0:
                    coarse_w = fine_grid_info['vertical_velocity'][::self.refinement_factor, ::self.refinement_factor]
                    coarse_p = fine_grid_info['pressure_perturbation'][::self.refinement_factor, ::self.refinement_factor]
                    self.vertical_velocity[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_w
                    self.pressure_perturbation[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_p
                if self.enable_baroclinic:
                    coarse_density = fine_grid_info['density'][::self.refinement_factor, ::self.refinement_factor]
                    self.density[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_density
                if self.enable_barotropic:
                    coarse_shear = fine_grid_info['velocity_shear'][::self.refinement_factor, ::self.refinement_factor]
                    self.velocity_shear[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_shear
                if self.enable_mli and grid_spacing_km < 10.0:
                    coarse_buoyancy = fine_grid_info['buoyancy'][::self.refinement_factor, ::self.refinement_factor]
                    self.buoyancy[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_buoyancy
                if self.enable_pv_frontogenesis:
                    coarse_pv = fine_grid_info['potential_vorticity'][::self.refinement_factor, ::self.refinement_factor]
                    self.potential_vorticity[j_start:j_start+patch_shape[0], i_start:i_start+patch_shape[1]] = coarse_pv
        
        self.update_visualization()

    def update_visualization(self):
        self.image.fill(Qt.white)
        painter = QPainter(self.image)
        if self.vorticity is not None:
            norm_vorticity = (self.vorticity - np.min(self.vorticity)) / (np.max(self.vorticity) - np.min(self.vorticity) + 1e-10)
            grad_x = (np.roll(self.vorticity, -1, axis=1) - np.roll(self.vorticity, 1, axis=1)) / (2 * 2/self.grid_size)
            grad_y = (np.roll(self.vorticity, -1, axis=0) - np.roll(self.vorticity, 1, axis=0)) / (2 * 2/self.grid_size)
            grad_magnitude = np.sqrt(grad_x**2 + grad_y**2 + 1e-10)
            norm_grad = (grad_magnitude - np.min(grad_magnitude)) / (np.max(grad_magnitude) - np.min(grad_magnitude) + 1e-10)
            for i in range(self.grid_size):
                for j in range(self.grid_size):
                    value = norm_vorticity[i, j]
                    grad = norm_grad[i, j]
                    r = int(255 * grad if grad > 0.7 else 0)
                    g = int(255 * (1 - value) * (1 - grad))
                    b = int(255 * value * (1 - grad))
                    if self.enable_nonhydrostatic and self.grid_spacing_km < 1.0:
                        w = self.vertical_velocity[i, j]
                        if w > 0:
                            b = int(b * (1 + 0.5 * w))
                        else:
                            r = int(r * (1 - 0.5 * w))
                    if self.enable_baroclinic:
                        rho = (self.density[i, j] - 1025.0) / 0.1
                        g = int(g * (1 + 0.2 * np.clip(rho, -1, 1)))
                    if self.enable_barotropic:
                        shear = self.velocity_shear[i, j] / 0.1
                        r = int(r * (1 + 0.2 * np.clip(shear, -1, 1)))
                    if self.enable_mli and self.grid_spacing_km < 10.0:
                        b_mli = self.buoyancy[i, j] / 0.1
                        b = int(b * (1 + 0.3 * np.clip(b_mli, -1, 1)))
                    if self.enable_pv_frontogenesis:
                        pv = self.potential_vorticity[i, j] / 1.0
                        r = int(r * (1 + 0.3 * np.clip(pv, -1, 1)))
                        b = int(b * (1 + 0.3 * np.clip(pv, -1, 1)))
                    color = QColor(max(0, min(255, r)), max(0, min(255, g)), max(0, min(255, b)))
                    painter.setPen(QPen(color))
                    painter.setBrush(QBrush(color))
                    x = int(400 * i / self.grid_size)
                    y = int(400 * j / self.grid_size)
                    painter.drawRect(x, y, 4, 4)
            
            # Visualize fine grids
            for fine_grid_info in self.fine_grids:
                fine_vorticity = fine_grid_info['vorticity']
                x0, y0 = fine_grid_info['center']
                radius = fine_grid_info['radius']
                norm_fine = (fine_vorticity - np.min(fine_vorticity)) / (np.max(fine_vorticity) - np.min(fine_vorticity) + 1e-10)
                grad_fine_x = (np.roll(fine_vorticity, -1, axis=1) - np.roll(fine_vorticity, 1, axis=1)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                grad_fine_y = (np.roll(fine_vorticity, -1, axis=0) - np.roll(fine_vorticity, 1, axis=0)) / (2 * 2/(self.fine_grid_size * self.refinement_factor))
                grad_fine = np.sqrt(grad_fine_x**2 + grad_fine_y**2 + 1e-10)
                norm_grad_fine = (grad_fine - np.min(grad_fine)) / (np.max(grad_fine) - np.min(grad_fine) + 1e-10)
                for i in range(self.fine_grid_size * self.refinement_factor):
                    for j in range(self.fine_grid_size * self.refinement_factor):
                        value = norm_fine[i, j]
                        grad = norm_grad_fine[i, j]
                        r = int(255 * grad if grad > 0.7 else 0)
                        g = int(255 * (1 - value) * (1 - grad))
                        b = int(255 * value * (1 - grad))
                        if self.enable_nonhydrostatic and self.grid_spacing_km < 1.0:
                            w = fine_grid_info['vertical_velocity'][i, j]
                            if w > 0:
                                b = int(b * (1 + 0.5 * w))
                            else:
                                r = int(r * (1 - 0.5 * w))
                        if self.enable_baroclinic:
                            rho = (fine_grid_info['density'][i, j] - 1025.0) / 0.1
                            g = int(g * (1 + 0.2 * np.clip(rho, -1, 1)))
                        if self.enable_barotropic:
                            shear = fine_grid_info['velocity_shear'][i, j] / 0.1
                            r = int(r * (1 + 0.2 * np.clip(shear, -1, 1)))
                        if self.enable_mli and self.grid_spacing_km < 10.0:
                            b_mli = fine_grid_info['buoyancy'][i, j] / 0.1
                            b = int(b * (1 + 0.3 * np.clip(b_mli, -1, 1)))
                        if self.enable_pv_frontogenesis:
                            pv = fine_grid_info['potential_vorticity'][i, j] / 1.0
                            r = int(r * (1 + 0.3 * np.clip(pv, -1, 1)))
                            b = int(b * (1 + 0.3 * np.clip(pv, -1, 1)))
                        color = QColor(max(0, min(255, r)), max(0, min(255, g)), max(0, min(255, b)))
                        painter.setPen(QPen(color))
                        painter.setBrush(QBrush(color))
                        x_fine = int(400 * (x0 - radius + (i / (self.fine_grid_size * self.refinement_factor) * 2 * radius) + 1) / 2)
                        y_fine = int(400 * (y0 - radius + (j / (self.fine_grid_size * self.refinement_factor) * 2 * radius) + 1) / 2)
                        painter.drawRect(x_fine, y_fine, 2, 2)
                
                painter.setPen(QPen(Qt.red, 2))
                painter.setBrush(Qt.NoBrush)
                x_start = int(400 * (x0 - radius + 1) / 2)
                y_start = int(400 * (y0 - radius + 1) / 2)
                width = int(400 * 2 * radius / 2)
                painter.drawRect(x_start, y_start, width, width)
            
            # Draw streamlines
            painter.setPen(QPen(Qt.black, 1))
            for _ in range(10):
                x0 = np.random.uniform(-1, 1)
                y0 = np.random.uniform(-1, 1)
                points = []
                for _ in range(20):
                    i = int((x0 + 1) * self.grid_size / 2)
                    j = int((y0 + 1) * self.grid_size / 2)
                    if 0 <= i < self.grid_size and 0 <= j < self.grid_size:
                        u, v = self.compute_induced_velocity(x0, y0, self.eddy_centers, self.eddy_strengths, 0.0, 0.0)
                        points.append((int(400 * (x0 + 1) / 2), int(400 * (y0 + 1) / 2)))
                        x0 += u * 0.01
                        y0 += v * 0.01
                    else:
                        break
                if len(points) > 1:
                    for k in range(len(points) - 1):
                        painter.drawLine(points[k][0], points[k][1], points[k + 1][0], points[k + 1][1])
        painter.end()
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.drawImage(0, 0, self.image)

class OceanicEddyAndFrontWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        try:
            self.setWindowTitle("Oceanic Eddies and Fronts Simulation")
            self.setGeometry(200, 200, 600, 1000)
            self.setFont(QFont("Consolas", 10))
            
            self.setStyleSheet("""
                QDialog, QWidget { 
                    background-color: #ECEFF1; 
                    color: #000000; 
                    border: 1px solid #B0BEC5;
                    border-radius: 10px;
                }
                QLineEdit { 
                    background-color: #FFFFFF; 
                    color: #000000; 
                    border: 2px solid #B0BEC5; 
                    border-radius: 5px;
                    padding: 5px;
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QPushButton { 
                    background-color: #008080; 
                    color: #FFFFFF; 
                    border: 2px solid #006666; 
                    border-radius: 8px;
                    padding: 8px; 
                    font-family: Consolas;
                    font-size: 10pt;
                    font-weight: bold;
                }
                QPushButton:hover { 
                    background-color: #006666; 
                }
                QPushButton:pressed { 
                    background-color: #004C4C; 
                    border: 2px inset #006666; 
                }
                QPushButton:disabled { 
                    background-color: #B0BEC5; 
                    color: #78909C; 
                }
                QLabel { 
                    color: #000000; 
                    padding: 5px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QGroupBox { 
                    background-color: #FFFFFF; 
                    border: 1px solid #B0BEC5; 
                    border-radius: 8px;
                    margin-top: 12px; 
                    padding: 10px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QGroupBox::title { 
                    subcontrol-origin: margin; 
                    subcontrol-position: top center; 
                    padding: 5px; 
                    background-color: #ECEFF1; 
                    color: #000000; 
                    font-family: Consolas;
                    font-size: 10pt;
                    font-weight: bold;
                }
                QCheckBox { 
                    color: #000000; 
                    padding: 5px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
            """)
            
            self.setup_ui()
            self.timer = QTimer()
            self.timer.timeout.connect(self.update_simulation)
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(False)
            self.export_button.setEnabled(False)
            logging.debug("OceanicEddyAndFrontWindow initialized")
        except Exception as e:
            logging.error(f"OceanicEddyAndFrontWindow initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize: {str(e)}")

    def setup_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)
        self.setLayout(main_layout)
        
        self.simulation_widget = EddySimulationWidget()
        self.simulation_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(self.simulation_widget, stretch=3)
        
        control_group = QGroupBox("Eddy Simulation Parameters")
        control_group.setStyleSheet("QGroupBox { background-color: #FFFFFF; box-shadow: 0px 2px 5px rgba(0,0,0,0.2); }")
        form_layout = QFormLayout()
        form_layout.setHorizontalSpacing(10)
        form_layout.setVerticalSpacing(8)
        
        self.inputs = {}
        params = [
            ("Grid Size", "100", "Coarse grid resolution (10–200)"),
            ("Grid Spacing (km)", "1.0", "Physical grid spacing (0.1–10.0 km)"),
            ("Eddy Strength", "1.0", "Vortex intensity (0.5–2.0)"),
            ("Eddy Radius", "0.2", "Vortex size (0.1–0.3)"),
            ("Number of Eddies", "3", "Number of vortices (1–10)"),
            ("Vorticity Diffusion", "0.01", "Diffusion rate (0.005–0.02)"),
            ("Rotation Rate (rad/s)", "0.1", "Spin speed (0.05–0.2)"),
            ("Background Flow (m/s)", "0.05", "Flow speed (0.02–0.1)"),
            ("Decay Rate (1/s)", "0.01", "Dissipation rate (0.005–0.02)"),
            ("Refinement Threshold", "0.5", "Vorticity threshold for refinement (0.1–1.0)"),
            ("Refinement Factor", "4", "Fine grid resolution factor (2–8)"),
            ("Strain Rate (1/s)", "0.02", "Strain for frontogenesis (0.01–0.05)"),
            ("Coriolis Parameter (1/s)", "1e-4", "Coriolis parameter (1e-5–1e-3)"),
            ("Rossby Number", "0.1", "Rossby number for ageostrophic effects (0.05–0.5)"),
            ("Baroclinic Shear (1/s)", "0.01", "Shear for baroclinic instability (0.005–0.02)"),
            ("Stratification Param", "0.1", "Stratification strength (0.05–0.2)"),
            ("Mixed-Layer Depth (m)", "50.0", "Mixed-layer depth for MLI (10–100 m)"),
            ("Buoyancy Gradient (m/s²)", "0.01", "Buoyancy gradient for MLI (0.005–0.02)"),
            ("PV Gradient Factor", "0.05", "PV gradient strength for frontogenesis (0.01–0.1)")
        ]
        
        for label, default, tooltip in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Consolas", 10))
            edit.setMinimumWidth(100)
            edit.setToolTip(tooltip)
            form_layout.addRow(QLabel(label), edit)
            self.inputs[label] = edit
        
        self.enable_nonhydrostatic = QCheckBox("Enable Non-Hydrostatic Dynamics")
        self.enable_nonhydrostatic.setFont(QFont("Consolas", 10))
        self.enable_nonhydrostatic.setToolTip("Enable non-hydrostatic effects for scales < 1 km")
        form_layout.addRow(self.enable_nonhydrostatic)
        
        self.enable_baroclinic = QCheckBox("Enable Baroclinic Instability")
        self.enable_baroclinic.setFont(QFont("Consolas", 10))
        self.enable_baroclinic.setToolTip("Enable baroclinic instability from density gradients")
        form_layout.addRow(self.enable_baroclinic)
        
        self.enable_barotropic = QCheckBox("Enable Barotropic Instability")
        self.enable_barotropic.setFont(QFont("Consolas", 10))
        self.enable_barotropic.setToolTip("Enable barotropic instability from velocity shear")
        form_layout.addRow(self.enable_barotropic)
        
        self.enable_mli = QCheckBox("Enable Mixed-Layer Instability")
        self.enable_mli.setFont(QFont("Consolas", 10))
        self.enable_mli.setToolTip("Enable submesoscale mixed-layer instability for scales < 10 km")
        form_layout.addRow(self.enable_mli)
        
        self.enable_pv_frontogenesis = QCheckBox("Enable PV Frontogenesis")
        self.enable_pv_frontogenesis.setFont(QFont("Consolas", 10))
        self.enable_pv_frontogenesis.setToolTip("Enable potential vorticity frontogenesis")
        form_layout.addRow(self.enable_pv_frontogenesis)
        
        control_group.setLayout(form_layout)
        main_layout.addWidget(control_group, stretch=1)
        
        button_layout = QHBoxLayout()
        button_layout.setSpacing(8)
        
        self.start_button = QPushButton("Start")
        self.start_button.setFont(QFont("Consolas", 10, QFont.Bold))
        self.start_button.setMinimumHeight(40)
        self.start_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.start_button.clicked.connect(self.start_simulation)
        button_layout.addWidget(self.start_button)
        
        self.pause_button = QPushButton("Pause")
        self.pause_button.setFont(QFont("Consolas", 10, QFont.Bold))
        self.pause_button.setMinimumHeight(40)
        self.pause_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.pause_button.clicked.connect(self.pause_simulation)
        button_layout.addWidget(self.pause_button)
        
        self.reset_button = QPushButton("Reset")
        self.reset_button.setFont(QFont("Consolas", 10, QFont.Bold))
        self.reset_button.setMinimumHeight(40)
        self.reset_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.reset_button.clicked.connect(self.reset_simulation)
        button_layout.addWidget(self.reset_button)
        
        self.initialize_button = QPushButton("Initialize")
        self.initialize_button.setFont(QFont("Consolas", 10, QFont.Bold))
        self.initialize_button.setMinimumHeight(40)
        self.initialize_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.initialize_button.clicked.connect(self.initialize_simulation)
        button_layout.addWidget(self.initialize_button)
        
        self.export_button = QPushButton("Export")
        self.export_button.setFont(QFont("Consolas", 10, QFont.Bold))
        self.export_button.setMinimumHeight(40)
        self.export_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.export_button.clicked.connect(self.export_simulation)
        button_layout.addWidget(self.export_button)
        
        self.close_button = QPushButton("Close")
        self.close_button.setFont(QFont("Consolas", 10, QFont.Bold))
        self.close_button.setMinimumHeight(40)
        self.close_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.close_button.clicked.connect(self.close)
        button_layout.addWidget(self.close_button)
        
        main_layout.addLayout(button_layout)

    def initialize_simulation(self):
        try:
            grid_size = int(self.inputs["Grid Size"].text())
            grid_spacing_km = float(self.inputs["Grid Spacing (km)"].text())
            eddy_strength = float(self.inputs["Eddy Strength"].text())
            eddy_radius = float(self.inputs["Eddy Radius"].text())
            num_eddies = int(self.inputs["Number of Eddies"].text())
            refinement_threshold = float(self.inputs["Refinement Threshold"].text())
            refinement_factor = int(self.inputs["Refinement Factor"].text())
            coriolis_param = float(self.inputs["Coriolis Parameter (1/s)"].text())
            rossby_number = float(self.inputs["Rossby Number"].text())
            baroclinic_shear = float(self.inputs["Baroclinic Shear (1/s)"].text())
            stratification_param = float(self.inputs["Stratification Param"].text())
            mixed_layer_depth = float(self.inputs["Mixed-Layer Depth (m)"].text())
            buoyancy_gradient = float(self.inputs["Buoyancy Gradient (m/s²)"].text())
            pv_gradient_factor = float(self.inputs["PV Gradient Factor"].text())
            enable_nonhydrostatic = self.enable_nonhydrostatic.isChecked()
            enable_baroclinic = self.enable_baroclinic.isChecked()
            enable_barotropic = self.enable_barotropic.isChecked()
            enable_mli = self.enable_mli.isChecked()
            enable_pv_frontogenesis = self.enable_pv_frontogenesis.isChecked()
            
            
            if grid_size < 10 or grid_size > 200:
                raise ValueError("Grid Size must be between 10 and 200")
            if grid_spacing_km <= 0 or grid_spacing_km > 10:
                raise ValueError("Grid Spacing must be between 0.1 and 10.0 km")
            if eddy_strength <= 0 or eddy_strength > 2.0:
                raise ValueError("Eddy Strength must be between 0.5 and 2.0")
            if eddy_radius <= 0 or eddy_radius > 0.3:
                raise ValueError("Eddy Radius must be between 0.1 and 0.3")
            if num_eddies < 1 or num_eddies > 10:
                raise ValueError("Number of Eddies must be between 1 and 10")
            if refinement_threshold <= 0 or refinement_threshold > 1.0:
                raise ValueError("Refinement Threshold must be between 0.1 and 1.0")
            if refinement_factor < 2 or refinement_factor > 8:
                raise ValueError("Refinement Factor must be between 2 and 8")
            if coriolis_param < 1e-5 or coriolis_param > 1e-3:
                raise ValueError("Coriolis Parameter must be between 1e-5 and 1e-3")
            if rossby_number < 0.05 or rossby_number > 0.5:
                raise ValueError("Rossby Number must be between 0.05 and 0.5")
            if baroclinic_shear < 0.005 or baroclinic_shear > 0.02:
                raise ValueError("Baroclinic Shear must be between 0.005 and 0.02")
            if stratification_param < 0.05 or stratification_param > 0.2:
                raise ValueError("Stratification Param must be between 0.05 and 0.2")
            if mixed_layer_depth < 10 or mixed_layer_depth > 100:
                raise ValueError("Mixed-Layer Depth must be between 10 and 100 m")
            if buoyancy_gradient < 0.005 or buoyancy_gradient > 0.02:
                raise ValueError("Buoyancy Gradient must be between 0.005 and 0.02 m/s²")
            if pv_gradient_factor < 0.01 or pv_gradient_factor > 0.1:
                raise ValueError("PV Gradient Factor must be between 0.01 and 0.1")
                
            self.simulation_widget.initialize_field(grid_size, eddy_strength, eddy_radius, num_eddies, 
                                                  refinement_threshold, refinement_factor, coriolis_param, 
                                                  rossby_number, enable_nonhydrostatic, grid_spacing_km, 
                                                  enable_baroclinic, enable_barotropic, baroclinic_shear, 
                                                  stratification_param, enable_mli, mixed_layer_depth, 
                                                  buoyancy_gradient, enable_pv_frontogenesis, pv_gradient_factor)
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.export_button.setEnabled(True)
            logging.debug("Eddy simulation initialized with ageostrophic, non-hydrostatic, baroclinic, barotropic, mixed-layer instability, and PV frontogenesis dynamics")
            
        except ValueError as e:
            logging.error(f"Initialization failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", f"Invalid input: {str(e)}")
        except Exception as e:
            logging.error(f"Unexpected initialization error: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Unexpected error: {str(e)}")

    def start_simulation(self):
        try:
            if self.simulation_widget.vorticity is None:
                self.initialize_simulation()
            vorticity_diffusion = float(self.inputs["Vorticity Diffusion"].text())
            rotation_rate = float(self.inputs["Rotation Rate (rad/s)"].text())
            background_flow = float(self.inputs["Background Flow (m/s)"].text())
            decay_rate = float(self.inputs["Decay Rate (1/s)"].text())
            refinement_threshold = float(self.inputs["Refinement Threshold"].text())
            strain_rate = float(self.inputs["Strain Rate (1/s)"].text())
            coriolis_param = float(self.inputs["Coriolis Parameter (1/s)"].text())
            rossby_number = float(self.inputs["Rossby Number"].text())
            enable_nonhydrostatic = self.enable_nonhydrostatic.isChecked()
            grid_spacing_km = float(self.inputs["Grid Spacing (km)"].text())
            enable_baroclinic = self.enable_baroclinic.isChecked()
            enable_barotropic = self.enable_barotropic.isChecked()
            baroclinic_shear = float(self.inputs["Baroclinic Shear (1/s)"].text())
            stratification_param = float(self.inputs["Stratification Param"].text())
            enable_mli = self.enable_mli.isChecked()
            mixed_layer_depth = float(self.inputs["Mixed-Layer Depth (m)"].text())
            buoyancy_gradient = float(self.inputs["Buoyancy Gradient (m/s²)"].text())
            enable_pv_frontogenesis = self.enable_pv_frontogenesis.isChecked()
            pv_gradient_factor = float(self.inputs["PV Gradient Factor"].text())
            
            self.timer.timeout.disconnect()
            self.timer.timeout.connect(lambda: self.simulation_widget.update_simulation(
                vorticity_diffusion, rotation_rate, background_flow, decay_rate, 
                refinement_threshold, strain_rate, coriolis_param, rossby_number, 
                enable_nonhydrostatic, grid_spacing_km, enable_baroclinic, enable_barotropic, 
                baroclinic_shear, stratification_param, enable_mli, mixed_layer_depth, 
                buoyancy_gradient, enable_pv_frontogenesis, pv_gradient_factor))
            self.timer.start(100)
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(True)
            self.reset_button.setEnabled(True)
            self.export_button.setEnabled(True)
            logging.debug("Eddy simulation started")
        except Exception as e:
            logging.error(f"Start simulation failed: {str(e)}")
            QMessageBox.critical(self, "Start Error", f"Failed to start simulation: {str(e)}")

    def pause_simulation(self):
        try:
            self.timer.stop()
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.export_button.setEnabled(True)
            logging.debug("Eddy simulation paused")
        except Exception as e:
            logging.error(f"Pause simulation failed: {str(e)}")
            QMessageBox.critical(self, "Pause Error", f"Failed to pause simulation: {str(e)}")

    def reset_simulation(self):
        try:
            self.timer.stop()
            self.simulation_widget.vorticity = None
            self.simulation_widget.vertical_velocity = None
            self.simulation_widget.pressure_perturbation = None
            self.simulation_widget.density = None
            self.simulation_widget.velocity_shear = None
            self.simulation_widget.buoyancy = None
            self.simulation_widget.potential_vorticity = None
            self.simulation_widget.time = 0.0
            self.simulation_widget.eddy_centers = []
            self.simulation_widget.eddy_strengths = []
            self.simulation_widget.eddy_radii = []
            self.simulation_widget.fine_grids = []
            self.simulation_widget.enable_nonhydrostatic = False
            self.simulation_widget.enable_baroclinic = False
            self.simulation_widget.enable_barotropic = False
            self.simulation_widget.enable_mli = False
            self.simulation_widget.enable_pv_frontogenesis = False
            self.simulation_widget.image.fill(Qt.white)
            self.simulation_widget.update()
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.export_button.setEnabled(False)
            logging.debug("Eddy simulation reset")
        except Exception as e:
            logging.error(f"Reset failed: {str(e)}")
            QMessageBox.critical(self, "Reset Error", f"Failed to reset simulation: {str(e)}")

    def export_simulation(self):
        try:
            if self.simulation_widget.vorticity is None:
                raise ValueError("No simulation data to export. Please initialize the simulation first.")
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            filename = f"oceanic_eddies_output_{timestamp}.txt"
            data = self.simulation_widget.get_simulation_data()
            data.insert(0, f"Oceanic Eddies and Fronts Simulation Output")
            data.insert(1, f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S %Z')}")
            data.insert(2, "-" * 50)
            with open(filename, 'w') as f:
                for line in data:
                    f.write(line + '\n')
            logging.debug(f"Simulation data exported to {filename}")
            QMessageBox.information(self, "Export Successful", f"Simulation data exported to {filename}")
        except ValueError as e:
            logging.error(f"Export failed: {str(e)}")
            QMessageBox.warning(self, "Export Error", f"Invalid input: {str(e)}")
        except Exception as e:
            logging.error(f"Unexpected export error: {str(e)}")
            QMessageBox.critical(self, "Export Error", f"Unexpected error during export: {str(e)}")
    def update_simulation(self):
        try:
            vorticity_diffusion = float(self.inputs["Vorticity Diffusion"].text())
            rotation_rate = float(self.inputs["Rotation Rate (rad/s)"].text())
            background_flow = float(self.inputs["Background Flow (m/s)"].text())
            decay_rate = float(self.inputs["Decay Rate (1/s)"].text())
            refinement_threshold = float(self.inputs["Refinement Threshold"].text())
            strain_rate = float(self.inputs["Strain Rate (1/s)"].text())
            coriolis_param = float(self.inputs["Coriolis Parameter (1/s)"].text())
            rossby_number = float(self.inputs["Rossby Number"].text())
            baroclinic_shear = float(self.inputs["Baroclinic Shear (1/s)"].text())
            stratification_param = float(self.inputs["Stratification Param"].text())
            mixed_layer_depth = float(self.inputs["Mixed-Layer Depth (m)"].text())
            buoyancy_gradient = float(self.inputs["Buoyancy Gradient (m/s²)"].text())
            enable_nonhydrostatic = self.enable_nonhydrostatic.isChecked()
            enable_baroclinic = self.enable_baroclinic.isChecked()
            enable_barotropic = self.enable_barotropic.isChecked()
            enable_mli = self.enable_mli.isChecked()
            enable_pv_frontogenesis = self.enable_pv_frontogenesis.isChecked()
            grid_spacing_km = float(self.inputs["Grid Spacing (km)"].text())
            
            if vorticity_diffusion < 0.005 or vorticity_diffusion > 0.02:
                raise ValueError("Vorticity Diffusion must be between 0.005 and 0.02")
            if rotation_rate < 0.05 or rotation_rate > 0.2:
                raise ValueError("Rotation Rate must be between 0.05 and 0.2")
            if background_flow < 0.02 or background_flow > 0.1:
                raise ValueError("Background Flow must be between 0.02 and 0.1")
            if decay_rate < 0.005 or decay_rate > 0.02:
                raise ValueError("Decay Rate must be between 0.005 and 0.02")
            if strain_rate < 0.01 or strain_rate > 0.05:
                raise ValueError("Strain Rate must be between 0.01 and 0.05")
            
            self.simulation_widget.update_simulation(vorticity_diffusion, rotation_rate, background_flow, 
                                                   decay_rate, refinement_threshold, strain_rate, 
                                                   coriolis_param, rossby_number, enable_nonhydrostatic, 
                                                   grid_spacing_km, enable_baroclinic, enable_barotropic, 
                                                   baroclinic_shear, stratification_param, enable_mli, 
                                                   mixed_layer_depth, buoyancy_gradient, enable_pv_frontogenesis)
        except ValueError as e:
            logging.error(f"Simulation update failed: {str(e)}")
            QMessageBox.warning(self, "Simulation Error", f"Invalid input: {str(e)}")
            self.pause_simulation()
        except Exception as e:
            logging.error(f"Unexpected simulation update error: {str(e)}")
            QMessageBox.critical(self, "Simulation Error", f"Unexpected error: {str(e)}")
            self.pause_simulation()

    def closeEvent(self, event):
        try:
            self.timer.stop()
            logging.debug("OceanicEddyAndFrontWindow closed")
            event.accept()
        except Exception as e:
            logging.error(f"Close event failed: {str(e)}")
            QMessageBox.critical(self, "Close Error", f"Failed to close window: {str(e)}")

