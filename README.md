# Software zur Modellierung und Simulation der Ozean-Atmosphäre-Kopplung (Software for modeling and simulating Ocean-Atmosphere Coupling)

_Die Software is geschrieben in Python-Programmiersprache und simuliert Wechselwirkungen zwischen Ozean und Atmosphäre mit einem Schwerpunkt auf physikalischen Prozessen wie Wärme-, Impuls-, Süßwasser- und CO₂-Flüssen, turbulenter Mischung, Grenzschichtdynamik, Wolkenmikrophysik, Wechselwirkungen zwischen Luft und Meer sowie ozeanischen Wirbeln und Fronten._

_The software, written in Python, simulates ocean-atmosphere interactions with a focus on physical processes such as heat, momentum, freshwater, and CO₂ fluxes, turbulent mixing, boundary layer dynamics, cloud microphysics, air-sea interactions, and oceanic eddies and fronts._

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
  
  $u_* = √(ρ_{air}·C_{d}·U^{2}/ρ_{water})$
  
  $z_{0} = α(u_{*}^{2} + 0.1(u_{ocean}^{2} + v_{ocean}^{2}))/g$

  $C_{d}' = C_{d}(1 + 0.1log10(z_{0}))$

  where $u_*$ is friction velocity, $C_{d}$ is drag coefficient, $U$ is wind speed, $ρ_{air} = 1.225 kg/m^{3}$, $ρ_{water} = 1025 kg/m^{3}$, $α = 0.018$, $g = 9.81 m/s^{2}$, and $z_{0}$ is roughness length.

- **Implementation**: `compute_sea_surface_roughness` clips $z_0$ to $[1e-6, 1e-2]$ and $C_{d}'$ to $[1e-4, 1e-2]$.

#### Momentum Flux
- **Purpose**: Computes wind stress on the ocean surface.
- **Equation**: $τ = ρ_{air}·C_{d}'·U^2$ where $C_{d}'$ is the adjusted drag coefficient.
- **Implementation**: `compute_momentum_flux` clips $τ$ to [-1e5, 1e5] $N/m^{2}$.

#### Heat Flux
- **Purpose**: Computes sensible and latent heat fluxes.
- **Equations**:

  $Q_{sensible} = ρ_{air}·C_{p}·k_{sensible}(T_{a} - T_{o})$

  $k_{sensible} = 0.01 * C_d'/C_d$

  $Q_{latent} = ρ_{air}·L_{v}·E·sign(T_{a} - T_{o})$

  $Q_{total} = Q_{sensible} + Q_{latent}$

  where $C_p = 1005 J/kg/K$, $L_v = 2.5e6 J/kg$, $E$ is evaporation rate, $T_a$ is atmosphere temperature, and $T_o$ is ocean temperature.

- **Implementation**: `compute_heat_flux` clips $T_a - T_o$ to $[-100, 100]$ $K$ and $Q_{total}$ to $[-1e6, 1e6]$ $W/m^{2}$.

#### Freshwater Flux
- **Purpose**: Computes net freshwater flux and salinity change, ensuring mass conservation.
- **Equations**:
  
  $P = P_{0}(1 + 0.5q/0.01)$

  $E = E_{0}·C_{freshwater}(1+0.1S/35)$

  $F = P - E$

  $dS/dt = -S·F/(ρ_{water}·H)$

  where $P_{0}$ is precipitation rate, $E_0$ is evaporation rate, $C_{freshwater}$ is the conservation coefficient, $S$ is salinity, $H$ = 1000 m is ocean depth, and $q$ is moisture.

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
- **Implementation**: `compute_co2_flux` clips outputs to $[-1e-3, 1e-3]$.

#### Moisture Advection
- **Purpose**: Computes moisture transport in the atmosphere.
- **Equation**: $Advection$ = $-u_{a}·∂q/∂x - v_{a}·∂q/∂y$ where $u_a$, $v_a$ are atmosphere velocities, and $∂q/∂x$, $∂q/∂y$ are computed via central differences.
- **Implementation**: `compute_moisture_advection` uses finite differences and clips output to $[-1e-4, 1e-4]$.

#### Turbulent Mixing
- **Purpose**: Models mixing driven by temperature gradients and wind.
- **Equation**:

  $M = k_{m}.u_*.(∂T/∂x + ∂T/∂y)$

  $u_* = sqrt(ρ_{air} * C_{d} * U^{2} / ρ_{water})$

  where $k_{m}$ is mixing coefficient, and gradients are computed via central differences.

- **Implementation**: `compute_turbulent_mixing` clips output to $[-1e3, 1e3]$.

#### Radiative Flux
- **Purpose**: Models solar and longwave radiation with greenhouse effects.
- **Equation**: $Q_{rad} = Q_{solar} - ε·σ·(T/300)^{4}·(1 + 0.1log(CO₂_{atm}/400))$ where $σ = 5.67e-8 W/m^{2}/K^{4}$, $ε$ is longwave coefficient, and $Q_{solar}$ is solar forcing.
- **Implementation**: `compute_radiative_flux` clips $T/300$ to $[0.8, 1.2]$, $log(CO₂_{atm}/400)$ to $[-1, 1]$, and output to $[-1e6, 1e6]$ $W/m^2$.

### 2. Ocean-Atmosphere Model 
The `OceanAtmosphereModel` integrates physical processes over time.

#### Initialization
- Sets up a 2D grid (N x N) with initial conditions:

  $T_o = T_{base} + ΔT·sin(2πy/N)$

  $T_a = T_o - 5$

  $S = S_0$

  $q = q_0$

  $CO₂_{ocean} = CO₂_{atm} = C_{0}$

  where $T_o$ is ocean temperature, $T_a$ is atmosphere temperature, $S$ is salinity, $q$ is moisture, and velocities are initialized to zero.

#### Time Stepping
- **Ocean Update**:
  
  $T_{o}^{n+1} = T_{o}^{n} + dt·(-u_{o}·∇T_{o} + k * ∇^{2}·T_{o} + (Q_{total} + Q_{rad})/(ρ_{water}·C_{p_{water}}·H))$

  $S^{n+1} = S^n + dt·dS/dt$

  $u_{o}^{n+1} = u_{o}^{n} + dt·(τ/(ρ_{water}H) - f·u_{o}^{⊥} + k·∇^{2}·u_{o})$

  where $k$ is diffusion coefficient, $f = 1e-4 s^{-1}$ is Coriolis parameter, and $∇^2$ is computed via central differences.

- **Atmosphere Update**:

  $T_{a}^{n+1} = T_{a}^{n} + dt/R·(-u_{a}·∇T_{a} + k·∇^{2}·T_{a} - (Q_{total} + Q_{rad})/(ρ_{air}·C_{p}^{air}·H_{atm}))$

  $q^{n+1} = q^n + dt/R * (Advection + k·∇^{2}·q + F/H_{atm})$

  $CO₂_{atm}^{n+1} = CO₂_{atm}^{n} + dt/R·F_{CO₂}^{atm}$

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
- **Euler Method**: $y^{n+1} = y^{n} + dt·f(y^{n})$ where $f$ is the flux or diffusivity rate.
- **RK4 Method**:

  $k_1 = f(y^n)$;

  $k_2 = f(y^n + dt/2 * k_1)$;

  $k_3 = f(y^n + dt/2 * k_2)$;

  $k_4 = f(y^n + dt * k_3)$;

  $y^{n+1} = y^n + dt/6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)$

- **Implementation**: `HighOrderTimeSteppingWindow.euler_step` and `rk4_step`.

### 4. Turbulent Mixing (`TurbulentMixing.py`)
- **Thermal Front**: $T_{front} = T_{cold} + (T_{hot} - T_{cold})·0.5·(1 + tanh((x-x_{0})/L))$ where $T_{cold} = 290 K$, $T_{hot} = 300 K$, $x_{0} = $N_{x}/2 + 10·sin(0.01 Ut)$, and $L = 5$.
- **Ekman Spiral**:

  $u = V_{0}·e^{z/d}·cos(π/4 + z/d)$

  $v = V_{0}·e^{z/d}·sin(π/4 + z/d)$

  $V_{0} = τ / (ρ·√(2νf))$

  $d = sqrt(2ν/f)$

  where $τ = 0.1U^{2}, $ν$ is mixing coefficient, and $f = 1e-4 s^{-1}$.

- **Vertical Mixing**:
  
  $dT/dt = K_{v}·d^{2}T/dz^{2}$

  $K_{v} = ν·e^{z/50}$
  
  where $d^{2}T/dz^{2}$ is computed via central differences.

- **Implementation**: `update_plot` blends front with ocean temperatures and computes Ekman velocities and temperature profiles.

### 5. Surface Layer Physics (`SurfaceLayerPhysics.py`)
- **Sensible Heat Flux**: $Q_s = ρ_{air}·C_p·C_h·U(T_o - T_a)$
- **Latent Heat Flux**: $Q_l = ρ_{air}·L_{v}·C_{e}·U·q$ where $q = 0.001$, $C_e$ is latent heat coefficient.
- **Wind Stress**: $τ = ρ_{air}·C_d·U^2$
- **Surface Currents**: $u_{s}^{n+1} = u_{s}^{n} + dt·(τ/1000)·N(0, ν)$ where $N$ is Gaussian noise scaled by mixing coefficient.
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
- **Potential Vorticity**: $q = (ζ + f)·(db/dz)/N^2$; $db/dz ≈ b/H_{m}$ where $ζ$ is vorticity, $f$ is Coriolis parameter, $b$ is buoyancy, $H_{m}$ is mixed layer depth, and $N^{2}$ is stratification parameter.
- **Frontogenesis**:

  $F = -|∇q|^{2}·D$

  $q^{n+1} = q^{n} + 0.1F$

  $dζ/dt = 0.1F$

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

  $∂φ/∂x ≈ (φ_{i+1,j} - φ_{i-1,j}) / (2·dx_{i,j})$

  - Uses central differences for spatial derivatives.
- **Diffusion**:
  $∂φ/∂t = k·∇^2·φ$

  $∇^{2}·φ ≈ (φ_{i+1,j} - 2φ_{i,j} + φ_{i-1,j})/dx_{i,j}^{2} + (φ_{i,j+1} - 2φ_{i,j} + φ_{i,j-1})/dy_{i,j}^{2}$

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

## Ozean-Atmosphären-Modell (Ocean-Atmosphere Model)

- **Initialization of Fields**:
   - Sets up a 2D grid ($NxN$) with initial conditions for ocean temperature ($T_o$), atmosphere temperature ($T_a$), salinity ($S$), ocean velocities ($u_{ocean}, v_{ocean}$), moisture ($q$), and CO₂ concentrations ($CO_{2}^{ocean}, CO_{2}^{atm}$).
   - Configures variable resolution grids, nested grids (optional), and adaptive mesh refinement (AMR).
   - Initializes physical and numerical parameters, including time scales for ocean and atmosphere.

- **Time Stepping**:
   - Advances the simulation using different time steps for ocean ($dt$·($ocean$ $time$ $scale$)) and atmosphere ($dt$·($atm$ $time$ $scale$)).
   - Updates fields by computing fluxes (heat, momentum, freshwater, CO₂), advection, diffusion, and turbulent mixing.
   - Applies numerical stability constraints (e.g., clipping temperatures to 250–350 K).

- **Physical Process Integration**:
   - Computes advection and diffusion for temperature, salinity, and CO₂ fields using finite difference methods.
   - Integrates fluxes and mixing from `TwoWayCoupling` to couple ocean and atmosphere dynamics.
   - Applies Coriolis effects (implicitly via momentum flux) and grid refinement based on temperature gradients or vorticity.

- **Grid Management**:
   - Uses `VariableResolutionGrid` for spatially variable resolution (finer near coasts).
   - Supports optional nested grids via `NestedGrid` for enhanced resolution in specific regions.
   - Employs `AdaptiveMeshRefinement` to dynamically refine the grid based on physical criteria.

- **Logging and Error Handling**:
   - Logs initialization and step progress using the `logging` module.
   - Includes exception handling to catch and log errors during simulation steps.

### _Logic_

- **Purpose**: Sets up the initial state of the ocean and atmosphere fields and configures numerical and physical parameters.
- **Process**:
  - Initializes 2D arrays for ocean temperature ($T_o$), atmosphere temperature ($T_a$), salinity ($S$), ocean velocities ($u_{ocean}$, $v_{ocean}$), moisture ($q$), and CO₂ concentrations ($CO_{2}^{ocean}$, $CO_{2}^{atm}$) with default values (e.g., $T_o = temp_{ocean}, T_a = temp_{atm}, S = 35.0 psu, q = 0.01, CO_{2}^{atm} = 400.0 ppm$).
  - Clips initial fields to physical ranges (e.g., $T_o, T_a$ to $[250, 350] K, S$ to $[30, 40] psu$).
  - Creates a `VariableResolutionGrid` instance to define spatially variable $dx$ and $dy$ (finer near coasts).
  - Optionally initializes a `NestedGrid` for high-resolution subdomains.
  - Sets up an `AdaptiveMeshRefinement` instance for dynamic grid refinement.
  - Configures a `TwoWayCoupling` instance with parameters for drag coefficient, wind speed, precipitation, evaporation, solar forcing, longwave coefficient, mixing coefficient, and CO₂ transfer coefficient.
  - Applies different time scales for ocean ($dt_{ocean} = dt$·($ocean$ $time$ $scale$)) and atmosphere ($dt_{atm} = dt$·($atm$ $time$ $scale$)), with defaults of 1.0 and 0.1, respectively.

### _Time Stepping (step method)_
- **Purpose**: Advances the simulation by one time step, updating all fields.
- **Process**:
  1. Copies current fields ($T_o, T_a, S, u_{ocean}, v_{ocean}, q, (CO_{2}^{ocean}, (CO_{2}^{atm}$) to avoid in-place modifications.
  2. Computes a refinement mask using `AdaptiveMeshRefinement.compute_refinement` based on temperature gradients or vorticity.
  3. Updates fields on the nested grid (if enabled) using `NestedGrid.update`.
  4. Defines time-varying atmosphere velocities:
     - $u_{atm} = velocity_{adv}·cos(2π·step·dt/$ $total$ $time$)
       ```
       u_atm = adv_velocity * cos(2 * π * step * dt / total_time)
       ```
     - $v_{atm} = velocity_{adv}·sin(2π·step·dt/$ $total$ $time$)
       ```
       v_atm = adv_velocity * sin(2 * π * step * dt / total_time)
       ```
  5. Computes fluxes and mixing using `TwoWayCoupling`:
     - Heat flux ($Q$) using $T_a, T_o, u_{ocean}, v_{ocean}$.
     - Radiative flux ($R_{ocean}, R_{atm}$) using $T_o$, $T_a$, and $CO_{2}^{atm}$.
     - Freshwater flux ($dS/dt$, $F_{freshwater}$) using $S$ and $q$.
     - Momentum flux ($τ$) using wind speed, $u_{ocean}$, $v_{ocean}$.
     - Moisture advection ($M_{adv}$) using $q, dx, dy, u_{atm}, v_{atm}$.
     - Turbulent mixing ($mix_{ocean}, mix_{atm}$) using $T_o, T_a, S, dx, dy$, wind speed.
     - CO₂ flux ($F_{CO_{2}}^{ocean}$, $F_{CO_{2}}^{atm}$) using $CO_{2}^{ocean}, CO_{2}^{atm}$.
  6. Updates ocean velocities:
     - $u_{new} = u_{ocean} + dt_{ocean}·τ/ρ_{water}$
       ```
       u_new = u_ocean + ocean_dt * τ / ρ_water
       ```
     - $v_{new} = v_{ocean} + dt_{ocean}·τ/ρ_{water}$
       ```
       v_new = v_ocean + ocean_dt * τ / ρ_water
       ```
  7. Computes advection for $T_o$, $T_a$, $S$, $CO_{2}^{ocean}$, $CO_{2}^{atm}$ using `compute_advection`.
  8. Computes diffusion for $T_o$, $T_a$, $S$ using `compute_diffusion`.
  9. Updates fields using a semi-implicit scheme ($α = 0.5$):
     - $(T_{o_{new}} = T_o + dt_{ocean}·(α·(Q/C_o + R_{ocean}/C_o - adv_{ocean} + diff_{ocean} + mix_{ocean}/C_o)+(1 - α)·(Q/C_o+R_{ocean}/C_o))$
       
       ```
       T_o_new = T_o + ocean_dt * (α * (Q / C_o + R_ocean / C_o - adv_ocean + diff_ocean + mix_ocean / C_o) + (1 - α) * (Q / C_o + R_ocean / C_o))
       ```
     - $T_{a_{new}} = T_a + dt_{atm}·(α·(-Q/C_a+R_{atm}/C_a-adv_{atm} + diff_{atm} + mix_{atm}/C_a) + (1 - α)·(-Q/C_a + R_{atm}/C_a))$
 
       ```
       T_a_new = T_a + atm_dt * (α * (-Q / C_a + R_atm / C_a - adv_atm + diff_atm + mix_atm / C_a) + (1 - α) * (-Q / C_a + R_atm / C_a))
       ```
     - $S_{new} = S + dt·(dS/dt + diff_{salinity} + mix_{ocean})$
    
       ```
       S_new = S + dt * (dS/dt + diff_salinity + mix_ocean)
       ```
     - $q_{new} = q + dt·(M_{adv} - F_{freshwater})$
    
       ```
       q_new = q + dt * (M_adv - F_freshwater)
       ```
     - $C_{o_{new}} = CO_{2}^{ocean} + dt·(F_{CO_{2}}^{ocean} + adv_{(CO_{2})_{o}})$
    
       ```
       C_o_new = CO₂_ocean + dt * (F_co2_ocean + adv_co2_o)
       ```
     - $C_{a_{new}} = CO_{2}^{atm} + dt·(F_{CO_{2}^{atm}} + adv_{(CO_{2})_{a}})$
    
       ```
       C_a_new = CO₂_atm + dt * (F_co2_atm + adv_co2_a)
       ```
  10. Clips updated fields to physical ranges (e.g., $T_o, T_a$ to $[250, 350]$ $K, u_ocean, v_ocean$ to $[-10, 10] m/s$, $q$ to $[0, 0.05]$, $CO_{2}^{atm}$ to $[200, 1000] ppm$).
  11. Applies AMR to refine $T_o$, $T_a$, and $S$ where the refinement mask is active.
  12. Returns the current time, updated $T_o$, $T_a$, and refinement mask.

### _Advection and Diffusion_
- **Advection (`compute_advection`)**:
  - Computes advection for a field ($T$) using velocities $u, v$:
    - $∂T/∂t = -u·(∂T/∂x) - v·(∂T/∂y)$
      
      ```
      ∂T/∂t = -u * ∂T/∂x - v * ∂T/∂y
      ```
    - $∂T/∂x ≈ (T_{i+1,j} - T_{i-1,j})/(2·dx_{i,j})$
   
      ```
      ∂T/∂x ≈ (T_(i+1,j) - T_(i-1,j)) / (2 * dx_i,j)
      ```
    - $∂T/∂y ≈ (T_{i,j+1} - T_{i,j-1})/(2·dy_{i,j})$
   
      ```
      ∂T/∂y ≈ (T_(i,j+1) - T_(i,j-1)) / (2 * dy_i,j)
      ```
  - Uses central differences and clips output to $[-1e5, 1e5]$.
- **Diffusion (`compute_diffusion`)**:
  - Computes diffusion for a field ($T$):
    - $∂T/∂t = D·∇^{2}·T$
    - $∇^{2}·T ≈ (T_{i+1,j} - 2·T_{i,j} + T_{i-1,j}) / dx_{i,j}^{2} + (T_{i,j+1} - 2·T_{i,j} + T_{i,j-1}) / dy_{i,j}^{2}$
   
    ```
    ∇²T ≈ (T_(i+1,j) - 2*T_i,j + T_(i-1,j)) / dx_i,j² + (T_(i,j+1) - 2*T_i,j + T_(i,j-1)) / dy_i,j²
    ```
    - $D = 1e-6$ (diffusion coefficient)
  - Clips output to $[-1e5, 1e5]$.

### _Ocean Temperature Update_
- Updates ocean temperature ($T_o$) based on heat flux, radiative flux, advection, diffusion, and turbulent mixing.
- $(T_{o})^{n+1} = (T_o)^{n} + dt_{ocean}·(α·(Q/C_o + R_{ocean} / C_o - adv_{ocean} + diff_{ocean} + mix_{ocean} /C_o) + (1 - α)·(Q/C_o + R_{ocean} / C_o))$

  ```
  T_o^(n+1) = T_o^n + ocean_dt * (α * (Q / C_o + R_ocean / C_o - adv_ocean + diff_ocean + mix_ocean / C_o) + (1 - α) * (Q / C_o + R_ocean / C_o))
  ```
  
  where:
  - $Q$ is heat flux from `TwoWayCoupling.compute_heat_flux` (W/m²).
  - $R_{ocean}$ is radiative flux from `TwoWayCoupling.compute_radiative_flux` (W/m²).
  - $adv_{ocean}$ is advection from `compute_advection` (K/s).
  - $diff_{ocean}$ is diffusion from `compute_diffusion` (K/s).
  - $mix_{ocean}$ is turbulent mixing from `TwoWayCoupling.compute_turbulent_mixing` (W/m²).
  - $C_o$ is ocean heat capacity (J/kg/K).
  - $α = 0.5$ is the semi-implicit factor.
  - $dt_{ocean} = dt$·($ocean$ $time$ $scale$).
- **Implementation**: Clips terms to $[-1e3, 1e3]$ and $T_o$ to $[250, 350] K$.

### _Atmosphere Temperature Update_
- Updates atmosphere temperature ($T_a$) based on heat flux, radiative flux, advection, diffusion, and turbulent mixing.
- $(T_a)^{n+1} = (T_a)^{n} + dt_{atm}·(α·(-Q/C_a + R_{atm}/C_a - adv_{atm} + diff_{atm} + mix_{atm} / C_a) + (1 - α)·(-Q / C_a + R_{atm} / C_a))$

  ```
  T_a^(n+1) = T_a^n + atm_dt * (α * (-Q / C_a + R_atm / C_a - adv_atm + diff_atm + mix_atm / C_a) + (1 - α) * (-Q / C_a + R_atm / C_a))
  ```
  where:
  - $Q$ is heat flux (negative for atmosphere due to coupling).
  - $R_{atm}$ is radiative flux.
  - $adv_{atm}$ is advection.
  - $diff_{atm}$ is diffusion.
  - $mix_{atm}$ is turbulent mixing.
  - $C_a$ is atmosphere heat capacity (J/kg/K).
  - $dt_{atm} = dt$·($atm$ $time$ $scale$).
- **Implementation**: Clips terms to $[-1e3, 1e3]$ and $T_a$ to $[250, 350] K$.

### _Salinity Update_
- Updates salinity ($S$) based on freshwater flux, diffusion, and turbulent mixing.
- $S^{n+1} = S^n + dt·(dS/dt + diff_{salinity} + mix_{ocean})$

  ```
  S^(n+1) = S^n + dt * (dS/dt + diff_salinity + mix_ocean)
  ```
  where:
  - $dS/dt$ is salinity change from `TwoWayCoupling.compute_freshwater_flux`.
  - $diff_{salinity}$ is diffusion.
  - $mix_{ocean}$ is turbulent mixing.
- **Implementation**: Clips terms to $[-1e-2, 1e-2]$ and $S$ to $[30, 40] psu$.

### _Moisture Update_
- Updates atmospheric moisture ($q$) based on advection and freshwater flux.
- $q^{n+1} = q^n + dt·(M_{adv} - F_{freshwater})$

  ```
  q^(n+1) = q^n + dt * (M_adv - F_freshwater)
  ```
  where:
  - $M_{adv}$ is moisture advection from `TwoWayCoupling.compute_moisture_advection`.
  - $F_{freshwater}$ is freshwater flux.
- **Implementation**: Clips terms to $[-1e-4, 1e-4]$ and $q$ to $[0, 0.05]$.

### _CO₂ Concentration Updates_
- Updates ocean and atmosphere CO₂ concentrations ($(CO_{2})^{ocean}$), ($CO_{2})^{atm}$) based on CO₂ flux and advection.
- $(CO_{2}^{ocean})^{n+1} = (CO_{2}^{ocean})^{n} + dt·((F_{CO_{2}})^{ocean} + (adv_{CO_{2}})_{o})$
  ```
  CO₂_ocean^(n+1) = CO₂_ocean^n + dt * (F_co2_ocean + adv_co2_o)
  ```
  $(CO_{{2}^{atm}})^{n+1} = (CO_{2}^{atm})^{n} + dt·((F_{CO_{2}})^{atm} + (adv_{CO_{2}})_{a})$
  where:
  - $(F_{CO_{2}})^{ocean}$, $(F_{CO_{2}})^{atm}$ are CO₂ fluxes from `TwoWayCoupling.compute_co2_flux`.
  - $(adv_{CO_{2}})ₒ)$, $(adv_{CO_{2}})_{a})$ are advection terms.
- **Implementation**: Clips terms to $[-1e-2, 1e-2]$, $CO_{2}^{ocean}$ to $[0, 10]$, and $CO_{{2}^{atm}}$ to $[200, 1000] ppm$.

### _Ocean Velocity Update_
- Updates ocean velocities ($u_{ocean}, v_{ocean}$) based on momentum flux.
- $u_{ocean}^{n+1} = u_{ocean}^{n} + dt_{ocean}·τ/ρ_{water}$

  ```
  u_ocean^(n+1) = u_ocean^n + ocean_dt * τ / ρ_water
  ```
  $v_{ocean}^{n+1} = v_{ocean}^{n} + dt_{ocean}·τ/ρ_{water}$
  ```
  v_ocean^(n+1) = v_ocean^n + ocean_dt * τ / ρ_water
  ```
  where:
  - $τ$ is momentum flux from `TwoWayCoupling.compute_momentum_flux` (N/m²).
  - $ρ_{water} = 1025 kg/m^{3}$ (from `TwoWayCoupling`).
- **Implementation**: Clips velocities to [-10, 10] m/s.

### _Advection_
- Computes advection for any field ($T, S, CO_{2}$).
- $∂T/∂t = -u·(∂T/∂x) - v·(∂T/∂y)$
  ```
  ∂T/∂t = -u * ∂T/∂x - v * ∂T/∂y
  ```
  $∂T/∂x ≈ (T_{i+1,j} - T_{i-1,j})/(2·dx_{i,j})$
  ```
  ∂T/∂x ≈ (T_(i+1,j) - T_(i-1,j)) / (2 * dx_i,j)
  ```
  $∂T/∂y ≈ (T_{i,j+1} - T_{i,j-1})/(2·dy_{i,j})$
  ```
  ∂T/∂y ≈ (T_(i,j+1) - T_(i,j-1)) / (2 * dy_i,j)
  ```
- **Implementation**: Uses central differences, supports array-based or scalar velocities, and clips output to [-1e5, 1e5].

### _Diffusion_
- Computes diffusion for temperature or salinity fields.
- $∂T/∂t = D·∇^{2}·T$

  $∇^{2}·T ≈ (T_{i+1,j} - 2T_{i,j} + T_{i-1,j})/dx_{i,j}^{2} + (T_{i,j+1} - 2T_{i,j} + T_{i,j-1})/dy_{i,j}^{2}$
  ```
  ∇²T ≈ (T_(i+1,j) - 2*T_i,j + T_(i-1,j))/(dx_i,j)² + (T_(i,j+1) - 2*T_i,j + T_(i,j-1))/(dy_i,j)²
  ```
  where $D = 1e-6 m^2/s$
- **Implementation**: Uses central differences and clips output to $[-1e5, 1e5]$.

---

![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(1).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(2).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(3).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(4).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(5).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(6).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(7).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(8).png)
![](https://github.com/KMORaza/Ocean-Atmospheric_Coupling_Software/blob/main/Ocean-Atmosphere%20Coupling%20Model%20Software/screenshot/screen%20(9).png)




















