# Software zur Modellierung und Simulation der Ozean-Atmosphäre-Kopplung (Software for modeling and simulating Ocean-Atmosphere Coupling)

The software, written in Python, simulates ocean-atmosphere interactions with a focus on physical processes such as heat, momentum, freshwater, and CO₂ fluxes, turbulent mixing, boundary layer dynamics, cloud microphysics, air-sea interactions, and oceanic eddies and fronts. 

## Funktionalitäten (Functionalities)

The ocean-atmosphere coupling aoftware is a numerical simulation software designed to model interactions between the ocean and atmosphere. It supports 2D simulations with variable resolution grids, adaptive mesh refinement, and two-way coupling for physical interactions. 

1. **Core Simulation**:
   - Simulates coupled ocean-atmosphere dynamics on a 2D grid with variable time steps for ocean and atmosphere.
   - Computes temperature, salinity, velocity, moisture, and CO₂ concentrations, accounting for advection, diffusion, and fluxes.

2. **Physical Processes**:
   - Models momentum, heat, freshwater, and CO₂ fluxes between ocean and atmosphere.
   - Simulates turbulent mixing, boundary layer schemes (Bulk and KPP), cloud microphysics, air-sea interactions, and surface layer physics.
   - Supports potential vorticity (PV) frontogenesis for oceanic eddies and fronts.

3. **Numerical Methods**:
   - Implements Euler and Runge-Kutta 4 (RK4) time-stepping methods for boundary layer schemes.
   - Uses finite difference methods for advection and diffusion.
   - Employs adaptive mesh refinement and variable resolution grids for computational efficiency.

4. **Visualization and Analysis**:
   - Visualizes ocean and atmosphere temperature fields, mean temperatures over time, and refinement regions.
   - Provides specialized analyses for turbulent mixing, surface layer physics, boundary layer schemes, air-sea interactions, and cloud microphysics.
   - Tracks thermal fronts, Ekman spirals, heat fluxes, wind stresses, and cloud properties.

5. **Parameter Control**:
   - Allows configuration of physical parameters (e.g., drag coefficient, wind speed, mixing coefficient) and numerical parameters (e.g., grid size, time step).
   - Supports conservation constraints for freshwater and CO₂ fluxes.

### Komponenten (Modules)
- **MainApp.py**: Orchestrates the simulation, integrating all modules.
- **Model.py**: Core simulation engine (`OceanAtmosphereModel`) managing state updates.
- **ControlPanel.py**: Manages simulation parameters.
- **TwoWayCoupling.py**: Computes ocean-atmosphere fluxes and turbulent mixing.
- **VariableResolution.py**: Defines variable-resolution grids for coastal regions.
- **AdaptiveMeshRefinement.py**: Implements dynamic grid refinement.
- **TurbulentMixing.py**: Analyzes turbulent mixing and Ekman spiral.
- **SurfaceLayerPhysics.py**: Models surface layer dynamics (heat fluxes, wind-driven currents).
- **HighOrderTimeStepping.py**: Analyzes high-order time-stepping for boundary layer schemes.
- **BoundaryLayerSchemes.py**: Implements Bulk and K-Profile Parameterization (KPP) schemes.
- **AirSeaInteraction.py**: Analyzes air-sea fluxes and coupling strength.
- **CloudMicroPhysics.py**: Simulates cloud formation and precipitation.
- **PotentialVorticityFrontogenesis.py**: Models PV and frontogenesis effects.
- **OceanicEddyAndFront.py**: Simulates oceanic eddies and fronts with non-hydrostatic, baroclinic, barotropic, and mixed-layer instabilities.
- **PlotWidget.py**: Visualizes simulation results (temperature, vorticity, refinement grids).
- **ConsoleWidget.py**: Logs simulation events.

## Simulationslogik (Simulation Logic)

1. **Initialization** (`Model.py`):
   - Sets up a 2D grid with initial conditions for ocean temperature, salinity, velocities, atmosphere temperature, moisture, and CO₂ concentrations.
   - Initializes a `VariableResolutionGrid` (`VariableResolution.py`) with finer resolution near coasts.
   - Configures a `TwoWayCoupling` instance (`TwoWayCoupling.py`) for flux calculations and an `AdaptiveMeshRefinement` instance (`AdaptiveMeshRefinement.py`) for grid refinement.

2. **Time Stepping** (`Model.py`):
   - Advances the simulation using different time steps for ocean (dt) and atmosphere (dt / time_scale_ratio).
   - Updates fields in the following order:
     - Computes heat, momentum, freshwater, and CO₂ fluxes using `TwoWayCoupling`.
     - Applies advection and diffusion for temperature, salinity, and moisture.
     - Updates velocities with momentum flux and Coriolis force.
     - Refines the grid using `AdaptiveMeshRefinement` based on temperature gradients or vorticity.
     - Applies numerical stability constraints (e.g., clipping temperatures to 250–350 K).

3. **Specialized Analyses**:
   - **Turbulent Mixing** (`TurbulentMixing.py`): Simulates a thermal front and vertical mixing, computing the Ekman spiral and temperature profiles.
   - **Surface Layer Physics** (`SurfaceLayerPhysics.py`): Models wind-driven surface currents and heat fluxes (sensible and latent).
   - **Boundary Layer Schemes** (`BoundaryLayerSchemes.py`): Analyzes Bulk and KPP schemes for heat flux or diffusivity.
   - **High-Order Time Stepping** (`HighOrderTimeStepping.py`): Compares Euler and RK4 methods for boundary layer schemes.
   - **Air-Sea Interaction** (`AirSeaInteraction.py`): Computes momentum and heat fluxes across the air-sea interface.
   - **Cloud Microphysics** (`CloudMicroPhysics.py`): Models cloud formation and precipitation based on moisture and temperature.
   - **Potential Vorticity Frontogenesis** (`PotentialVorticityFrontogenesis.py`): Simulates PV dynamics, likely used by `OceanicEddyAndFront.py` for eddy and front analyses.

4. **Visualization** (`PlotWidget.py`):
   - Displays ocean and atmosphere temperature fields as heatmaps, with contours for refinement regions and rectangles for nested grids.
   - Plots mean ocean and atmosphere temperatures over time as a time series.

5. **Parameter Input** (`ControlPanel.py`):
   - Provides parameters for initial conditions, physical constants, time scales, grid settings, and fluxes.
   - Triggers simulation steps and opens analysis windows for specific processes.

## Physikalische und Mathematische Modelle (Physics and Mathematical Models)

### 1. Two-Way Coupling 
The `TwoWayCoupling` class computes fluxes and mixing between the ocean and atmosphere, ensuring conservation where applicable.

#### Sea Surface Roughness
- **Purpose**: Adjusts drag coefficient based on wind and ocean currents.
- **Equation**:
  
  $u_* = sqrt(ρ_{air}C_{d}U^{2}$/$ρ_{water})$
  
  $z_{0} = α(u_{*}^{2} + 0.1(u_{ocean}^{2} + v_{ocean}^{2}))/g$

  $C_{d}' = C_{d}(1 + 0.1log10(z_{0}))$

  where $u_*$ is friction velocity, $C_{d}$ is drag coefficient, $U$ is wind speed, $ρ_{air} = 1.225 kg/m^{3}$, $ρ_{water} = 1025 kg/m^{3}$, $α = 0.018$, $g = 9.81 m/s^{2}, and $z_{0}$ is roughness length.

- **Implementation**: `compute_sea_surface_roughness` clips $z_0$ to [1e-6, 1e-2] and $C_{d}'$ to [1e-4, 1e-2].

#### Momentum Flux
- **Purpose**: Computes wind stress on the ocean surface.
- **Equation**: $τ = ρ_{air}*C_{d}'*U^2$ where $C_{d}'$ is the adjusted drag coefficient.
- **Implementation**: `compute_momentum_flux` clips $τ$ to [-1e5, 1e5] $N/m^{2}$.

#### Heat Flux
- **Purpose**: Computes sensible and latent heat fluxes.
- **Equations**:

  $Q_{sensible} = ρ_{air}C_{p}k_{sensible}(T_{a} - T_{o})$

  $k_{sensible} = 0.01 * C_d' / C_d$

  $Q_{latent} = ρ_{air}L_{v}E*sign(T_{a} - T_{o})$

  $Q_{total} = Q_{sensible} + Q_{latent}$

  where $C_p = 1005 J/kg/K$, $L_v = 2.5e6 J/kg$, $E$ is evaporation rate, $T_a$ is atmosphere temperature, and $T_o$ is ocean temperature.

- **Implementation**: `compute_heat_flux` clips $T_a - T_o$ to $[-100, 100]$ $K$ and $Q_{total}$ to $[-1e6, 1e6]$ $W/m^{2}$.

#### Freshwater Flux
- **Purpose**: Computes net freshwater flux and salinity change, ensuring mass conservation.
- **Equations**:
  
  $P = P_{0}*(1 + 0.5q/0.01)$

  $E = E_{0}*C_{freshwater}(1+0.1S/35)$

  $F = P - E$

  $dS/dt = -S*F/(ρ_{water}*H)$

  where $P_{0}$ is precipitation rate, $E_0$ is evaporation rate, $C_{freshwater} is the conservation coefficient, $S$ is salinity, $H$ = 1000 m is ocean depth, and $q$ is moisture.

- **Implementation**: `compute_freshwater_flux` clips $q/0.01$ to $[0, 2]$, $S/35$ to $[0.8, 1.2]$, and outputs to $[-1e-3, 1e-3]$.

#### CO₂ Flux
- **Purpose**: Computes CO₂ exchange between ocean and atmosphere, ensuring mass conservation.
- **Equations**:
  
  $pCO₂_{ocean} = CO₂_{ocean}/α$

  $pCO₂_{atm} = CO₂_{atm}$

  $F_{CO₂} = k_{CO₂}C_{CO₂}(pCO₂_{ocean}-pCO₂_{atm})$

  $F_{CO₂}^{ocean} = F_{CO₂}/H$

  $F_{CO₂}^{atm} = -F_{CO₂}/H_{atm}$

  where $α = 0.03$ is CO₂ solubility, $k_{CO₂}$ is transfer coefficient, $C_{CO₂}$ is conservation coefficient, $H = 1000 m$, and $H_{atm} = 10000 m$.
- **Implementation**: `compute_co2_flux` clips outputs to [-1e-3, 1e-3].

#### Moisture Advection
- **Purpose**: Computes moisture transport in the atmosphere.
- **Equation**: $Advection$ = $-u_a * ∂q/∂x - v_a * ∂q/∂y$ where $u_a$, $v_a$ are atmosphere velocities, and $∂q/∂x$, $∂q/∂y$ are computed via central differences.
- **Implementation**: `compute_moisture_advection` uses finite differences and clips output to [-1e-4, 1e-4].

#### Turbulent Mixing
- **Purpose**: Models mixing driven by temperature gradients and wind.
- **Equation**:

  $M = k_{m} * u_* (∂T/∂x + ∂T/∂y)$

  $u_* = sqrt(ρ_{air} * C_{d} * U^{2} / ρ_{water})$

  where $k_{m}$ is mixing coefficient, and gradients are computed via central differences.

- **Implementation**: `compute_turbulent_mixing` clips output to $[-1e3, 1e3]$.

#### Radiative Flux
- **Purpose**: Models solar and longwave radiation with greenhouse effects.
- **Equation**: $Q_{rad} = Q_{solar} - ε * σ * (T / 300)^4 * (1 + 0.1 * log(CO₂_{atm} / 400))$ where $σ = 5.67e-8 W/m^{2}/K^{4}$, $ε$ is longwave coefficient, and $Q_{solar}$ is solar forcing.
- **Implementation**: `compute_radiative_flux` clips $T/300$ to $[0.8, 1.2]$, $log(CO₂_{atm}/400)$ to $[-1, 1]$, and output to [-1e6, 1e6] $W/m^2$.

### 2. Ocean-Atmosphere Model 
The `OceanAtmosphereModel` integrates physical processes over time.

#### Initialization
- Sets up a 2D grid (N x N) with initial conditions:

  $T_o = T_{base} + ΔT*sin(2πy/N)$

  $T_a = T_o - 5$

  $S = S_0$

  $q = q_0$

  $CO₂_{ocean} = CO₂_{atm} = C_{0}$

  where $T_o$ is ocean temperature, $T_a$ is atmosphere temperature, $S$ is salinity, $q$ is moisture, and velocities are initialized to zero.

#### Time Stepping
- **Ocean Update**:
  
  $T_{o}^{n+1} = T_{o}^{n} + dt*(-u_{o}·∇T_{o} + k * ∇^{2}·T_{o} + (Q_{total} + Q_{rad}) / (ρ_{water}*C_{p_{water}} * H))$

  $S^{n+1} = S^n + dt*dS/dt$

  $u_{o}^{n+1} = u_{o}^{n} + dt*(τ/(ρ_{water}H) - f*u_{o}^{⊥} + k·∇^{2}·u_{o})$

  where $k$ is diffusion coefficient, $f = 1e-4 s^{-1}$ is Coriolis parameter, and $∇^2$ is computed via central differences.

- **Atmosphere Update**:

  $T_{a}^{n+1} = T_{a}^{n} + dt/R*(-u_{a}·∇T_{a} + k·∇^{2}·T_{a} - (Q_{total} + Q_{rad}) / (ρ_{air} * C_{p}^{air} * H_{atm}))$

  $q^{n+1} = q^n + dt/R * (Advection + k·∇^{2}·q + F/H_{atm})$

  $CO₂_{atm}^{n+1} = CO₂_{atm}^{n} + dt/R * F_{CO₂}^{atm}$

  where $R$ is the time scale ratio.

- **Numerical Stability**: Clips $T_o,$ $T_a$ to $[250, 350] K$, velocities to $[-0.5, 0.5] m/s$, and other fields to prevent divergence.

### 3. Boundary Layer Schemes
These modules model surface boundary layer processes using Bulk and KPP schemes.

#### Bulk Scheme
- **Heat Flux**: $Q = C_{h}U(T_o - T_a)$ where $C_h$ is sensible heat coefficient.
- **Implementation**: `BoundaryLayerSchemesWindow.compute_bulk_flux` and `HighOrderTimeSteppingWindow.compute_bulk`.

#### KPP Scheme
- **Diffusivity**: $K = k_{kpp}*(1-z/h)^{2}$ where $k_{kpp}$ is mixing coefficient, $z$ is depth, and $h$ is boundary layer depth.
- **Implementation**: `BoundaryLayerSchemesWindow.compute_kpp_diffusivity` and `HighOrderTimeSteppingWindow.compute_kpp`.

#### Time Stepping (`HighOrderTimeStepping.py`)
- **Euler Method**: $y^{n+1} = y^{n} + dt*f(y^{n})$ where $f$ is the flux or diffusivity rate.
- **RK4 Method**:

  $k_1 = f(y^n)$;

  $k_2 = f(y^n + dt/2 * k_1)$;

  $k_3 = f(y^n + dt/2 * k_2)$;

  $k_4 = f(y^n + dt * k_3)$;

  $y^{n+1} = y^n + dt/6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)$

- **Implementation**: `HighOrderTimeSteppingWindow.euler_step` and `rk4_step`.

### 4. Turbulent Mixing (`TurbulentMixing.py`)
- **Thermal Front**: $T_{front} = T_{cold} + (T_{hot} - T_{cold}) * 0.5 * (1 + tanh((x-x_{0})/L))$ where $T_{cold} = 290 K$, $T_{hot} = 300 K$, $x_{0} = $N_{x}/2 + 10 * sin(0.01 Ut)$, and $L = 5$.
- **Ekman Spiral**:

  $u = V_{0}·e^{z/d}·cos(π/4 + z/d)$

  $v = V_{0}·e^{z/d}·sin(π/4 + z/d)$

  $V_{0} = τ / (ρ*sqrt(2νf))$

  $d = sqrt(2ν/f)$

  where $τ = 0.1U^{2}, $ν$ is mixing coefficient, and $f = 1e-4 s^{-1}$.

- **Vertical Mixing**:
  
  $dT/dt = K_{v}·d^{2}T/dz^{2}$

  $K_{v} = ν·e^{z/50}$
  
  where $d^{2}T/dz^{2}$ is computed via central differences.

- **Implementation**: `update_plot` blends front with ocean temperatures and computes Ekman velocities and temperature profiles.

### 5. Surface Layer Physics (`SurfaceLayerPhysics.py`)
- **Sensible Heat Flux**: $Q_s = ρ_{air}·C_p·C_h·U(T_o - T_a)$
- **Latent Heat Flux**: $Q_l = ρ_{air} * L_v * C_e * U * q$ where $q = 0.001$, $C_e$ is latent heat coefficient.
- **Wind Stress**: $τ = ρ_{air}·C_d·U^2$
- **Surface Currents**: $u_{s}^{n+1} = u_{s}^{n} + dt*(τ/1000)*N(0, ν)$ where $N$ is Gaussian noise scaled by mixing coefficient.
- **Implementation**: `update_plot` computes fluxes and updates currents, clipped to $[-0.5, 0.5] m/s$.

### 6. Air-Sea Interaction (`AirSeaInteraction.py`)
- **Momentum Flux**: $τ = ρ_{air}·C_{d}·U^{2}$
- **Heat Flux**: $Q = ρ_{air}·C_{p}·C_{h}·U·(T_o - T_a)$
- **Implementation**: `update_plot` computes fluxes using `TwoWayCoupling` and visualizes as heatmaps.

### 7. Cloud Microphysics (`CloudMicroPhysics.py`)
- **Cloud Formation**:

  $C = 1$ if $q > q_{sat}(T_{a})$, else 0

  $q_{sat}(T) = 0.01·e^{0.06(T-273.15)}$

- **Precipitation**: $P = k_{p}·(q - q_{sat}(T_{a}))$ if $q > q_{sat}(T_{a})$, else 0 where $k_{p} = 0.01$.
- **Implementation**: `update_plot` computes cloud cover and precipitation based on moisture and temperature.

### 8. Potential Vorticity Frontogenesis 
- **Potential Vorticity**: $q = (ζ + f)*(db/dz)/N^2$; $db/dz ≈ b/H_{m}$ where $ζ$ is vorticity, $f$ is Coriolis parameter, $b$ is buoyancy, $H_{m}$ is mixed layer depth, and $N^{2}$ is stratification parameter.
- **Frontogenesis**:

  $F = -|∇q|^{2}·D$

  $q^{n+1} = q^{n} + 0.1*F$

  $dζ/dt = 0.1*F$

  where $D$ is deformation rate, and gradients are computed via central differences.
- **Implementation**: `compute_pv_frontogenesis` updates PV and returns vorticity tendency, clipped to [-0.1, 0.1].

### 9. Adaptive Mesh Refinement (`AdaptiveMeshRefinement.py`)
- **Refinement Criterion**: $|∇T|^2 = (∂T/∂x)^2 + (∂T/∂y)^2 > threshold$ or $ζ = ∂v/∂x - ∂u/∂y > vorticity_{threshold}$
- **Nested Grid**:
  - Places a finer grid (size $N/4$) at the region with maximum gradient or vorticity.
- **Implementation**: `refine_grid` computes gradients/vorticity and sets a refinement mask.

### 10. Variable Resolution Grid (`VariableResolution.py`)
- **Spatial Steps**: $dx_{i,j} = dy_{i,j} = 1 / coast factor$ if $i < N/4$ or $j < N/4$, else 1
- **Implementation**: `VariableResolutionGrid` initializes finer resolution near coasts.

## Algorithmen (Algorithms)

### Main Simulation Loop 
1. Initialize `OceanAtmosphereModel` with parameters from `ControlPanel`.
2. Start QTimer to call `step_simulation` every 100 ms.
3. In each step:
   - Call `OceanAtmosphereModel.step` to update fields.
   - Update `PlotWidget` with new temperatures, refinement mask, and nested grid parameters.
   - Log events to `ConsoleWidget`.

### Advection and Diffusion (`Model.py`)
- **Advection**:

  $∂φ/∂t = -u·∇φ$

  $∂φ/∂x ≈ (φ_{i+1,j} - φ_{i-1,j}) / (2*dx_{i,j})$

  - Uses central differences for spatial derivatives.
- **Diffusion**:
  $∂φ/∂t = k·∇^2·φ$

  $∇^{2}·φ ≈ (φ_{i+1,j} - 2φ_{i,j} + φ_{i-1,j})/dx_{i,j^{2}} + (φ_{i,j+1} - 2φ_{i,j} + φ_{i,j-1})/dy_{i,j^{2}}$

   - Applies diffusion to stabilize numerical solutions.

### Time Stepping for Boundary Layers (`HighOrderTimeStepping.py`)
- **Euler**: Single-step update using current rate.
- **RK4**: Computes four intermediate slopes ($k_1$, $k_2$, $k_3$, $k_4$) and combines them for higher accuracy.

### Turbulent Mixing Simulation (`TurbulentMixing.py`)
1. Initialize a thermal front and vertical grid ($z$) for the Ekman spiral.
2. Update temperature field: $T = 0.8T_{o} + 0.2T_{front}$.
3. Compute mixing field using `TwoWayCoupling.compute_turbulent_mixing`.
4. Update vertical temperature profile with diffusion.
5. Compute Ekman spiral velocities ($u, v$).
6. Visualize thermal front (contour) and vertical profiles (line plots).

### Surface Layer Physics Simulation (`SurfaceLayerPhysics.py`)
1. Compute sensible and latent heat fluxes using input parameters.
2. Compute wind stress and update surface currents with random mixing.
3. Visualize surface currents (heatmap) and time series of mean fluxes.

### Air-Sea Interaction Simulation (`AirSeaInteraction.py`)
1. Compute momentum and heat fluxes using `TwoWayCoupling`.
2. Visualize fluxes as heatmaps.

### Cloud Microphysics Simulation (`CloudMicroPhysics.py`)
1. Compute cloud cover based on moisture saturation.
2. Calculate precipitation rate for supersaturated regions.
3. Visualize cloud cover (heatmap) and precipitation (time series).

### Potential Vorticity Frontogenesis (`PotentialVorticityFrontogenesis.py`)
1. Initialize PV field using vorticity and buoyancy.
2. Compute PV gradients and frontogenesis term.
3. Update PV and return vorticity tendency.


