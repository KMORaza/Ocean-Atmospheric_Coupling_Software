import numpy as np

class PVFrontogenesis:
    def __init__(self, grid_size, grid_spacing_km, mixed_layer_depth, coriolis_param, stratification_param):
        """Initialize potential vorticity field."""
        self.grid_size = grid_size
        self.grid_spacing_km = grid_spacing_km
        self.mixed_layer_depth = mixed_layer_depth
        self.coriolis_param = coriolis_param
        self.stratification_param = stratification_param
        self.potential_vorticity = None

    def initialize_pv(self, vorticity, buoyancy):
        """Initialize potential vorticity based on vorticity and buoyancy."""
        if vorticity is None or buoyancy is None:
            self.potential_vorticity = np.zeros((self.grid_size, self.grid_size))
            return
        
        # Approximate vertical buoyancy gradient: db/dz ≈ buoyancy / mixed_layer_depth
        db_dz = buoyancy / self.mixed_layer_depth
        # PV = (vorticity + f) * (db/dz) / N^2
        self.potential_vorticity = (vorticity + self.coriolis_param) * db_dz / self.stratification_param
        self.potential_vorticity = np.clip(self.potential_vorticity, -1.0, 1.0)  # Prevent numerical instability

    def compute_pv_frontogenesis(self, vorticity, buoyancy, deformation_rate, dx):
        """Compute PV frontogenesis and its effect on vorticity."""
        if self.potential_vorticity is None:
            return np.zeros_like(vorticity)
        
        # Compute PV gradients
        grad_pv_x = (np.roll(self.potential_vorticity, -1, axis=1) - np.roll(self.potential_vorticity, 1, axis=1)) / (2 * dx)
        grad_pv_y = (np.roll(self.potential_vorticity, -1, axis=0) - np.roll(self.potential_vorticity, 1, axis=0)) / (2 * dx)
        
        # Deformation-driven frontogenesis: F = -|∇q|^2 * deformation_rate
        grad_pv_magnitude = np.sqrt(grad_pv_x**2 + grad_pv_y**2 + 1e-10)
        frontogenesis_term = -deformation_rate * grad_pv_magnitude**2
        
        # Update PV
        self.potential_vorticity += frontogenesis_term * 0.1  # Scale for stability
        self.potential_vorticity = np.clip(self.potential_vorticity, -1.0, 1.0)
        
        # Vorticity tendency from PV frontogenesis
        vorticity_tendency = 0.1 * frontogenesis_term  # Simplified coupling
        return np.clip(vorticity_tendency, -0.1, 0.1)
