import numpy as np
import logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QMessageBox, QGroupBox, QComboBox, QSlider
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QTimer, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
import matplotlib

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class AirSeaInteractionWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        logging.debug("Initializing AirSeaInteractionWindow")
        try:
            self.setWindowTitle("Enhanced Air-Sea Interaction Analysis")
            self.setGeometry(150, 150, 1000, 700)
            
            # Set Consolas font and 90s-style stylesheet
            self.setFont(QFont("Consolas", 10))
            self.setStyleSheet("""
                QMainWindow, QWidget { 
                    background-color: #C0C0C0; 
                    color: #000000; 
                    border: 2px solid #808080; 
                }
                QLineEdit { 
                    background-color: #FFFFFF; 
                    color: #000000; 
                    border: 2px inset #808080; 
                    padding: 2px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QPushButton { 
                    background-color: #008080; 
                    color: #FFFFFF; 
                    border: 2px outset #808080; 
                    padding: 4px; 
                    font-family: Consolas;
                    font-size: 10pt;
                    font-weight: bold;
                }
                QPushButton:hover { 
                    background-color: #006666; 
                }
                QPushButton:pressed { 
                    background-color: #004C4C; 
                    border: 2px inset #808080; 
                }
                QPushButton:disabled { 
                    background-color: #669999; 
                    color: #A0A0A0; 
                }
                QComboBox { 
                    background-color: #FFFFFF; 
                    color: #000000; 
                    border: 2px inset #808080; 
                    padding: 2px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QComboBox QAbstractItemView { 
                    background-color: #FFFFFF; 
                    color: #000000; 
                    border: 2px outset #808080; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QLabel { 
                    color: #000000; 
                    padding: 2px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QGroupBox { 
                    background-color: #C0C0C0; 
                    border: 2px groove #808080; 
                    margin-top: 10px; 
                    padding: 10px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QGroupBox::title { 
                    subcontrol-origin: margin; 
                    subcontrol-position: top left; 
                    padding: 0 3px; 
                    background-color: #C0C0C0; 
                    color: #000000; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QSlider { 
                    background-color: #C0C0C0; 
                    border: 2px inset #808080; 
                }
                QSlider::groove:horizontal { 
                    border: 1px solid #808080; 
                    height: 8px; 
                    background: #FFFFFF; 
                    margin: 2px 0; 
                }
                QSlider::handle:horizontal { 
                    background: #008080; 
                    border: 1px solid #808080; 
                    width: 18px; 
                    margin: -2px 0; 
                }
            """)
            
            # Initialize internal model parameters
            self.nx, self.ny = 100, 100  # Grid size
            self.dt = 0.1  # Time step (s)
            self.total_time = 100.0  # Total simulation time (s)
            self.dx = np.ones((self.ny, self.nx)) * 1000.0  # Spatial step (m)
            self.dy = np.ones((self.ny, self.nx)) * 1000.0  # Spatial step (m)
            self.rho_air = 1.225  # Air density (kg/m³)
            self.rho_water = 1025.0  # Water density (kg/m³)
            self.Cp_air = 1005.0  # Specific heat of air (J/kg/K)
            self.Lv = 2.5e6  # Latent heat of vaporization (J/kg)
            self.R = 287.0  # Gas constant for air (J/kg/K)
            self.Hco2 = 3.4e-2  # Henry's law constant for CO2 (mol/m³/atm)
            
            # Initialize synthetic data fields
            x = np.linspace(0, self.nx * 1000, self.nx)
            y = np.linspace(0, self.ny * 1000, self.ny)
            X, Y = np.meshgrid(x, y)
            self.ocean_temps = 290 + 10 * np.sin(2 * np.pi * X / (self.nx * 1000))  # Synthetic ocean temps (K)
            self.atm_temps = 295 + 5 * np.cos(2 * np.pi * Y / (self.ny * 1000))  # Synthetic atm temps (K)
            self.salinity = np.full((self.ny, self.nx), 35.0)  # Constant salinity (psu)
            self.moisture = np.full((self.ny, self.nx), 0.01)  # Constant moisture (kg/kg)
            self.u_ocean = 0.1 * np.cos(2 * np.pi * X / (self.nx * 1000))  # Synthetic ocean u-velocity (m/s)
            self.v_ocean = 0.1 * np.sin(2 * np.pi * Y / (self.ny * 1000))  # Synthetic ocean v-velocity (m/s)
            self.co2_ocean = np.full((self.ny, self.nx), 400.0)  # Ocean CO2 concentration (ppm)
            self.co2_atm = np.full((self.ny, self.nx), 410.0)  # Atmospheric CO2 concentration (ppm)
            
            # Default coupling parameters
            self.drag_coeff = 0.001
            self.wind_speed = 10.0
            self.wind_angle = 0.0  # Wind direction (degrees)
            self.ocean_current_speed = 0.1  # Ocean current speed (m/s)
            self.evap_rate = 1e-6
            self.precip_rate = 1e-6
            self.co2_transfer_coeff = 1e-5  # CO2 transfer coefficient (m/s)
            self.plot_x = self.nx // 2  # Default x-coordinate for time series
            self.plot_y = self.ny // 2  # Default y-coordinate for time series
            
            # Main widget and layout
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_layout = QHBoxLayout()
            main_widget.setLayout(main_layout)
            
            # Left side: Control panel
            left_layout = QVBoxLayout()
            control_group = QGroupBox("Air-Sea Interaction Parameters")
            form_layout = QFormLayout()
            form_layout.setHorizontalSpacing(10)
            
            self.inputs = {}
            params = [
                ("Drag Coefficient", str(self.drag_coeff)),
                ("Wind Speed (m/s)", str(self.wind_speed)),
                ("Wind Angle (degrees)", str(self.wind_angle)),
                ("Ocean Current Speed (m/s)", str(self.ocean_current_speed)),
                ("Evaporation Rate (kg/m²/s)", str(self.evap_rate)),
                ("Precipitation Rate (kg/m²/s)", str(self.precip_rate)),
                ("CO2 Transfer Coeff (m/s)", str(self.co2_transfer_coeff)),
                ("Plot X Coordinate", str(self.plot_x)),
                ("Plot Y Coordinate", str(self.plot_y))
            ]
            
            for label, default in params:
                edit = QLineEdit(default)
                edit.setFont(QFont("Consolas", 10))
                form_layout.addRow(QLabel(label), edit)
                self.inputs[label] = edit
            
            # Variable selection
            self.variable_selector = QComboBox()
            self.variable_selector.addItems(["Wind Stress", "Heat Flux", "Freshwater Flux", "CO2 Flux"])
            self.variable_selector.setFont(QFont("Consolas", 10))
            form_layout.addRow(QLabel("Variable"), self.variable_selector)
            
            # Animation speed slider
            self.speed_label = QLabel("Animation Speed (ms): 100")
            self.speed_label.setFont(QFont("Consolas", 10))
            self.speed_slider = QSlider(Qt.Horizontal)
            self.speed_slider.setMinimum(50)
            self.speed_slider.setMaximum(500)
            self.speed_slider.setValue(100)
            self.speed_slider.setTickInterval(50)
            self.speed_slider.setTickPosition(QSlider.TicksBelow)
            self.speed_slider.valueChanged.connect(self.update_speed)
            form_layout.addRow(self.speed_label, self.speed_slider)
            
            # Control buttons
            button_layout = QHBoxLayout()
            self.start_button = QPushButton("Start Simulation")
            self.start_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.start_button.clicked.connect(self.start_simulation)
            
            self.pause_button = QPushButton("Pause")
            self.pause_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.pause_button.clicked.connect(self.pause_simulation)
            self.pause_button.setEnabled(False)
            
            self.resume_button = QPushButton("Resume")
            self.resume_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.resume_button.clicked.connect(self.resume_simulation)
            self.resume_button.setEnabled(False)
            
            self.stop_button = QPushButton("Stop Simulation")
            self.stop_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.stop_button.clicked.connect(self.stop_simulation)
            self.stop_button.setEnabled(False)
            
            button_layout.addWidget(self.start_button)
            button_layout.addWidget(self.pause_button)
            button_layout.addWidget(self.resume_button)
            button_layout.addWidget(self.stop_button)
            
            form_layout.addRow(button_layout)
            control_group.setLayout(form_layout)
            left_layout.addWidget(control_group)
            left_layout.addStretch()
            
            # Right side: Plot widget
            self.figure = Figure(figsize=(8, 6), facecolor="#C0C0C0")
            self.canvas = FigureCanvas(self.figure)
            matplotlib.rc("axes", edgecolor="black", labelcolor="black", facecolor="#FFFFFF")
            matplotlib.rc("xtick", color="black")
            matplotlib.rc("ytick", color="black")
            matplotlib.rc("text", color="black")
            matplotlib.rc("font", family="Consolas", size=10)
            
            # Add layouts to main layout
            main_layout.addLayout(left_layout, 1)
            main_layout.addWidget(self.canvas, 3)
            
            # Initialize plot elements
            self.figure.clear()
            self.ax1 = self.figure.add_subplot(221)  # Hovmöller
            self.ax2 = self.figure.add_subplot(222)  # Heatmap
            self.ax3 = self.figure.add_subplot(212)  # Time series
            self.current_step = 0
            self.anim = None
            self.is_paused = False
            
            # Simulation setup
            self.time_steps = []
            self.data_history = []
            self.wind_stress = np.zeros((self.ny, self.nx))
            self.heat_flux = np.zeros((self.ny, self.nx))
            self.freshwater_flux = np.zeros((self.ny, self.nx))
            self.co2_flux = np.zeros((self.ny, self.nx))
            self.time_series_data = {
                "Wind Stress": [],
                "Heat Flux": [],
                "Freshwater Flux": [],
                "CO2 Flux": []
            }
            
            # Initialize timer
            self.timer = QTimer(self)
            self.timer.timeout.connect(self.start_animation)
            
            # Initialize plot
            self.update_plot(0)
            
            logging.debug("AirSeaInteractionWindow initialization complete")
        except Exception as e:
            logging.error(f"AirSeaInteractionWindow initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize air-sea interaction window: {str(e)}")
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
    
    def compute_momentum_flux(self, wind_speed, wind_angle, u_ocean, v_ocean):
        logging.debug("Computing momentum flux")
        try:
            drag_coeff = self.compute_sea_surface_roughness(wind_speed, u_ocean, v_ocean)
            wind_u = wind_speed * np.cos(np.radians(wind_angle))
            wind_v = wind_speed * np.sin(np.radians(wind_angle))
            rel_u = wind_u - u_ocean
            rel_v = wind_v - v_ocean
            speed_rel = np.sqrt(rel_u**2 + rel_v**2)
            tau = self.rho_air * drag_coeff * np.clip(speed_rel**2, 0, 1e3)
            tau_u = tau * rel_u / np.clip(speed_rel, 1e-6, None)
            tau_v = tau * rel_v / np.clip(speed_rel, 1e-6, None)
            return np.clip(tau, -1e5, 1e5), np.clip(tau_u, -1e5, 1e5), np.clip(tau_v, -1e5, 1e5)
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
            precip = self.precip_rate * (1 + 0.5 * np.clip(moisture / 0.01, 0, 2))
            evap = self.evap_rate * (1 + 0.1 * np.clip(salinity / 35.0, 0.8, 1.2))
            F = precip - evap
            dS_dt = -salinity * F / (self.rho_water * ocean_depth)
            return np.clip(dS_dt, -1e-3, 1e-3), np.clip(F, -1e-3, 1e-3)
        except Exception as e:
            logging.error(f"Freshwater flux computation failed: {str(e)}")
            raise
    
    def compute_co2_flux(self, co2_ocean, co2_atm, wind_speed, u_ocean, v_ocean):
        logging.debug("Computing CO2 flux")
        try:
            k_co2 = self.co2_transfer_coeff * (1 + 0.1 * np.clip(wind_speed / 10.0, 0.5, 2.0))
            delta_co2 = co2_ocean - self.Hco2 * co2_atm
            F_co2 = k_co2 * delta_co2
            return np.clip(F_co2, -1e-3, 1e-3)
        except Exception as e:
            logging.error(f"CO2 flux computation failed: {str(e)}")
            raise
    
    def update_speed(self, value):
        self.speed_label.setText(f"Animation Speed (ms): {value}")
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
            if not self.is_paused:
                self.anim = FuncAnimation(self.figure, self.update_plot, interval=value, 
                                       blit=False, cache_frame_data=False)
                self.canvas.draw()
        logging.debug(f"Animation speed updated to {value} ms")
    
    def start_simulation(self):
        logging.debug("Starting air-sea interaction simulation")
        try:
            self.drag_coeff = float(self.inputs["Drag Coefficient"].text())
            self.wind_speed = float(self.inputs["Wind Speed (m/s)"].text())
            self.wind_angle = float(self.inputs["Wind Angle (degrees)"].text())
            self.ocean_current_speed = float(self.inputs["Ocean Current Speed (m/s)"].text())
            self.evap_rate = float(self.inputs["Evaporation Rate (kg/m²/s)"].text())
            self.precip_rate = float(self.inputs["Precipitation Rate (kg/m²/s)"].text())
            self.co2_transfer_coeff = float(self.inputs["CO2 Transfer Coeff (m/s)"].text())
            self.plot_x = int(self.inputs["Plot X Coordinate"].text())
            self.plot_y = int(self.inputs["Plot Y Coordinate"].text())
            
            if (self.drag_coeff <= 0 or self.wind_speed < 0 or self.ocean_current_speed < 0 or 
                self.evap_rate < 0 or self.precip_rate < 0 or self.co2_transfer_coeff < 0):
                raise ValueError("Parameters must be non-negative; drag coefficient and CO2 transfer coefficient must be non-zero.")
            if not (0 <= self.plot_x < self.nx and 0 <= self.plot_y < self.ny):
                raise ValueError(f"Plot coordinates must be within grid (0-{self.nx-1}, 0-{self.ny-1}).")
            
            # Reset simulation data
            self.current_step = 0
            self.time_steps = []
            self.data_history = []
            self.wind_stress = np.zeros((self.ny, self.nx))
            self.heat_flux = np.zeros((self.ny, self.nx))
            self.freshwater_flux = np.zeros((self.ny, self.nx))
            self.co2_flux = np.zeros((self.ny, self.nx))
            self.time_series_data = {
                "Wind Stress": [],
                "Heat Flux": [],
                "Freshwater Flux": [],
                "CO2 Flux": []
            }
            
            # Update synthetic fields
            x = np.linspace(0, self.nx * 1000, self.nx)
            y = np.linspace(0, self.ny * 1000, self.ny)
            X, Y = np.meshgrid(x, y)
            self.ocean_temps = 290 + 10 * np.sin(2 * np.pi * X / (self.nx * 1000))
            self.atm_temps = 295 + 5 * np.cos(2 * np.pi * Y / (self.ny * 1000))
            self.u_ocean = self.ocean_current_speed * np.cos(2 * np.pi * X / (self.nx * 1000))
            self.v_ocean = self.ocean_current_speed * np.sin(2 * np.pi * Y / (self.ny * 1000))
            self.co2_ocean = 400 + 10 * np.sin(2 * np.pi * X / (self.nx * 1000))
            self.co2_atm = 410 + 5 * np.cos(2 * np.pi * Y / (self.ny * 1000))
            
            self.timer.start(self.speed_slider.value())
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(True)
            self.resume_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            self.is_paused = False
            logging.debug("Air-sea interaction simulation started")
        except ValueError as e:
            logging.error(f"Start air-sea interaction failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            logging.error(f"Start air-sea interaction failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to start air-sea interaction: {str(e)}")
    
    def pause_simulation(self):
        logging.debug("Pausing air-sea interaction simulation")
        try:
            if self.anim is not None:
                self.anim.event_source.stop()
                self.is_paused = True
            self.timer.stop()
            self.pause_button.setEnabled(False)
            self.resume_button.setEnabled(True)
            self.stop_button.setEnabled(True)
            logging.debug("Air-sea interaction simulation paused")
        except Exception as e:
            logging.error(f"Pause air-sea interaction failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to pause air-sea interaction: {str(e)}")
    
    def resume_simulation(self):
        logging.debug("Resuming air-sea interaction simulation")
        try:
            if self.is_paused:
                self.anim = FuncAnimation(self.figure, self.update_plot, interval=self.speed_slider.value(), 
                                       blit=False, cache_frame_data=False)
                self.timer.start(self.speed_slider.value())
                self.is_paused = False
            self.pause_button.setEnabled(True)
            self.resume_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            logging.debug("Air-sea interaction simulation resumed")
        except Exception as e:
            logging.error(f"Resume air-sea interaction failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to resume air-sea interaction: {str(e)}")
    
    def stop_simulation(self):
        logging.debug("Stopping air-sea interaction simulation")
        try:
            self.timer.stop()
            if self.anim is not None:
                self.anim.event_source.stop()
                self.anim = None
            self.is_paused = False
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.resume_button.setEnabled(False)
            self.stop_button.setEnabled(False)
            logging.debug("Air-sea interaction simulation stopped")
        except Exception as e:
            logging.error(f"Stop air-sea interaction failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to stop air-sea interaction: {str(e)}")
    
    def start_animation(self):
        if self.anim is None and not self.is_paused:
            self.anim = FuncAnimation(self.figure, self.update_plot, interval=self.speed_slider.value(), 
                                   blit=False, cache_frame_data=False)
            self.canvas.draw()
    
    def update_plot(self, frame):
        logging.debug(f"Updating air-sea interaction plot at step {self.current_step}")
        try:
            if self.current_step >= int(self.total_time / self.dt):
                self.stop_simulation()
                return
            
            variable = self.variable_selector.currentText()
            
            # Update synthetic fields with time-varying component
            x = np.linspace(0, self.nx * 1000, self.nx)
            y = np.linspace(0, self.ny * 1000, self.ny)
            X, Y = np.meshgrid(x, y)
            self.ocean_temps = 290 + 10 * np.sin(2 * np.pi * X / (self.nx * 1000) + 0.1 * self.current_step * self.dt)
            self.atm_temps = 295 + 5 * np.cos(2 * np.pi * Y / (self.ny * 1000) + 0.05 * self.current_step * self.dt)
            self.u_ocean = self.ocean_current_speed * np.cos(2 * np.pi * X / (self.nx * 1000) + 0.05 * self.current_step * self.dt)
            self.v_ocean = self.ocean_current_speed * np.sin(2 * np.pi * Y / (self.ny * 1000) + 0.05 * self.current_step * self.dt)
            self.co2_ocean = 400 + 10 * np.sin(2 * np.pi * X / (self.nx * 1000) + 0.05 * self.current_step * self.dt)
            self.co2_atm = 410 + 5 * np.cos(2 * np.pi * Y / (self.ny * 1000) + 0.05 * self.current_step * self.dt)
            
            # Compute fluxes
            wind_stress, tau_u, tau_v = self.compute_momentum_flux(self.wind_speed, self.wind_angle, self.u_ocean, self.v_ocean)
            self.wind_stress = wind_stress
            self.heat_flux = self.compute_heat_flux(self.atm_temps, self.ocean_temps, self.u_ocean, self.v_ocean)
            dS_dt, freshwater = self.compute_freshwater_flux(self.salinity, self.moisture)
            self.freshwater_flux = freshwater
            self.co2_flux = self.compute_co2_flux(self.co2_ocean, self.co2_atm, self.wind_speed, self.u_ocean, self.v_ocean)
            
            # Select data based on variable
            if variable == "Wind Stress":
                data = self.wind_stress
                unit = "N/m²"
                cmap = "viridis"
            elif variable == "Heat Flux":
                data = self.heat_flux
                unit = "W/m²"
                cmap = "coolwarm"
            elif variable == "CO2 Flux":
                data = self.co2_flux
                unit = "mol/m²/s"
                cmap = "magma"
            else:  # Freshwater Flux
                data = self.freshwater_flux
                unit = "kg/m²/s"
                cmap = "plasma"
            
            # Store data for Hovmöller and time series
            y_mid = self.ny // 2
            self.time_steps.append(self.current_step * self.dt)
            self.data_history.append(data[y_mid, :])
            self.time_series_data["Wind Stress"].append(self.wind_stress[self.plot_y, self.plot_x])
            self.time_series_data["Heat Flux"].append(self.heat_flux[self.plot_y, self.plot_x])
            self.time_series_data["Freshwater Flux"].append(self.freshwater_flux[self.plot_y, self.plot_x])
            self.time_series_data["CO2 Flux"].append(self.co2_flux[self.plot_y, self.plot_x])
            
            # Update plots
            self.ax1.clear()
            self.ax2.clear()
            self.ax3.clear()
            
            # Hovmöller diagram
            hov_data = np.array(self.data_history)
            extent = [0, self.nx, 0, len(self.time_steps) * self.dt]
            if len(self.data_history) >= 1:
                plot_data = hov_data if len(self.data_history) > 1 else hov_data.reshape(1, -1)
                im1 = self.ax1.imshow(plot_data, cmap=cmap, origin="lower", aspect="auto", extent=extent)
                self.ax1.set_title(f"{variable} Hovmöller (y={y_mid})")
                self.ax1.set_xlabel("X")
                self.ax1.set_ylabel("Time (s)")
                if not hasattr(self.ax1, 'colorbar') or not self.ax1.colorbar:
                    self.figure.colorbar(im1, ax=self.ax1, label=unit)
                    self.ax1.colorbar = True
            
            # Heatmap with contour lines and vector field
            x = np.arange(self.nx)
            y = np.arange(self.ny)
            X, Y = np.meshgrid(x, y)
            im2 = self.ax2.imshow(data, cmap=cmap, origin="lower")
            self.ax2.contour(X, Y, data, colors='black', linewidths=0.5)
            self.ax2.set_title(f"{variable} at t={self.current_step * self.dt:.2f}s")
            self.ax2.set_xlabel("X")
            self.ax2.set_ylabel("Y")
            if not hasattr(self.ax2, 'colorbar') or not self.ax2.colorbar:
                self.figure.colorbar(im2, ax=self.ax2, label=unit)
                self.ax2.colorbar = True
            
            # Wind stress vector field
            skip = max(1, self.nx // 10)
            X_sub, Y_sub = X[::skip, ::skip], Y[::skip, ::skip]
            U = tau_u[::skip, ::skip] / np.clip(wind_stress[::skip, ::skip], 1e-6, None)
            V = tau_v[::skip, ::skip] / np.clip(wind_stress[::skip, ::skip], 1e-6, None)
            self.ax2.quiver(X_sub, Y_sub, U, V, color='white', scale=0.5, width=0.002)
            
            # Mark selected point
            self.ax2.plot(self.plot_x, self.plot_y, 'r+', markersize=10)
            
            # Time series plot
            for var, values in self.time_series_data.items():
                if values:
                    self.ax3.plot(self.time_steps, values, label=var)
            self.ax3.set_title(f"Time Series at ({self.plot_x}, {self.plot_y})")
            self.ax3.set_xlabel("Time (s)")
            self.ax3.set_ylabel("Flux Values")
            self.ax3.legend()
            self.ax3.grid(True)
            
            self.figure.subplots_adjust(left=0.1, right=0.95, wspace=0.3, top=0.95, bottom=0.1, hspace=0.3)
            
            self.current_step += 1
            logging.debug(f"Air-sea interaction plot updated at step {self.current_step}")
            
            return [im1, im2]
        except ValueError as e:
            self.stop_simulation()
            logging.error(f"Air-sea interaction plot update failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            self.stop_simulation()
            logging.error(f"Air-sea interaction plot update failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to update air-sea interaction plot: {str(e)}")
    
    def closeEvent(self, event):
        logging.debug("Closing air-sea interaction window")
        try:
            self.timer.stop()
            if self.anim is not None:
                self.anim.event_source.stop()
                self.anim = None
            self.is_paused = False
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.resume_button.setEnabled(False)
            self.stop_button.setEnabled(False)
            event.accept()
        except Exception as e:
            logging.error(f"Close air-sea interaction window failed: {str(e)}")
            event.accept()