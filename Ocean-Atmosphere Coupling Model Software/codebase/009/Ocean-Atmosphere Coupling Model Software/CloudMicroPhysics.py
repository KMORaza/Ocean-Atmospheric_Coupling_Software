import numpy as np
import logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QMessageBox, QGroupBox, QComboBox, QSlider, QTextEdit, QScrollArea
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QTimer, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
import matplotlib
from datetime import datetime

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class CloudMicroPhysicsWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        logging.debug("Initializing CloudMicroPhysicsWindow")
        try:
            self.setWindowTitle("Enhanced Cloud Microphysics Analysis")
            self.setGeometry(200, 200, 1000, 700)
            
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
                QTextEdit { 
                    background-color: #FFFFFF; 
                    color: #000000; 
                    border: 2px inset #808080; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QScrollArea { 
                    background-color: #C0C0C0; 
                    border: 2px inset #808080; 
                }
            """)
            
            # Initialize internal model parameters
            self.n_droplets = 100  # Number of droplets to simulate
            self.dt = 0.1  # Time step (s)
            self.total_time = 100.0  # Total simulation time (s)
            self.rho_air = 1.225  # Air density (kg/m³)
            self.rho_water = 1000  # Water density (kg/m³)
            self.R = 287.0  # Gas constant for air (J/kg/K)
            self.Lv = 2.5e6  # Latent heat of vaporization (J/kg)
            self.Rv = 461.5  # Gas constant for water vapor (J/kg/K)
            self.c_p = 1004  # Specific heat capacity of air (J/kg/K)
            
            # Initialize synthetic data fields
            self.droplet_radii = np.random.normal(10.0, 1.0, self.n_droplets)  # Initial radii (μm)
            self.droplet_concentrations = np.full(self.n_droplets, 100.0)  # Initial concentration (cm⁻³)
            self.lwc = (4/3) * np.pi * (self.droplet_radii * 1e-6)**3 * 1e6 * self.droplet_concentrations  # g/m³
            self.supersaturation = 0.005  # Initial supersaturation (fraction)
            self.z = np.zeros(self.n_droplets)  # Height (m), initialized at cloud base
            self.rain_rate = 0.0  # Rain rate (mm/hr)
            
            # Default parameters
            self.updraft_vel = 1.0
            self.wind_speed = 10.0
            self.wind_angle = 0.0
            self.drag_coeff = 0.001
            self.air_temp = 293.15  # Air temperature (K)
            self.autoconversion_threshold = 20.0  # Droplet radius threshold for autoconversion (μm)
            self.wind_shear = 0.01  # Wind shear (s⁻¹)
            self.aerosol_concentration = 100.0  # Aerosol concentration (cm⁻³)
            self.aerosol_soluble_fraction = 0.5  # Soluble fraction for Köhler theory
            
            # Main widget and layout
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_layout = QHBoxLayout()
            main_widget.setLayout(main_layout)
            
            # Left side: Control panel and console
            left_layout = QVBoxLayout()
            control_group = QGroupBox("Cloud Microphysics Parameters")
            form_layout = QFormLayout()
            form_layout.setHorizontalSpacing(10)
            
            self.inputs = {}
            params = [
                ("Updraft Velocity (m/s)", str(self.updraft_vel)),
                ("Initial Supersaturation (%)", "0.5"),
                ("Initial Droplet Radius (μm)", "10.0"),
                ("Initial Droplet Concentration (cm⁻³)", "100.0"),
                ("Aerosol Concentration (cm⁻³)", str(self.aerosol_concentration)),
                ("Aerosol Soluble Fraction", str(self.aerosol_soluble_fraction)),
                ("Air Temperature (K)", str(self.air_temp)),
                ("Wind Speed (m/s)", str(self.wind_speed)),
                ("Wind Angle (degrees)", str(self.wind_angle)),
                ("Autoconversion Threshold (μm)", str(self.autoconversion_threshold)),
                ("Wind Shear (s⁻¹)", str(self.wind_shear))
            ]
            
            for label, default in params:
                edit = QLineEdit(default)
                edit.setFont(QFont("Consolas", 10))
                form_layout.addRow(QLabel(label), edit)
                self.inputs[label] = edit
            
            # Visualization option
            self.plot_selector = QComboBox()
            self.plot_selector.addItems(["Color by LWC", "Color by Droplet Size", "Droplet Size Histogram", "Vertical Profile"])
            self.plot_selector.setFont(QFont("Consolas", 10))
            form_layout.addRow(QLabel("Plot Style"), self.plot_selector)
            
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
            
            # Status label
            self.status_label = QLabel("Status: Ready")
            self.status_label.setFont(QFont("Consolas", 10))
            form_layout.addRow(QLabel("Simulation Status"), self.status_label)
            
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
            
            # Console for textual output
            console_group = QGroupBox("Simulation Output")
            console_layout = QVBoxLayout()
            self.console = QTextEdit()
            self.console.setReadOnly(True)
            self.console.setFont(QFont("Consolas", 10))
            console_layout.addWidget(self.console)
            console_group.setLayout(console_layout)
            
            # Wrap console in a scroll area
            console_scroll = QScrollArea()
            console_scroll.setWidgetResizable(True)
            console_scroll.setWidget(console_group)
            console_scroll.setMinimumHeight(150)
            console_scroll.setMaximumHeight(300)
            
            left_layout.addWidget(console_scroll)
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
            self.ax1 = self.figure.add_subplot(221)  # Scatterplot or histogram
            self.ax2 = self.figure.add_subplot(222)  # Histogram or vertical profile
            self.ax3 = self.figure.add_subplot(212)  # Time series
            self.current_step = 0
            self.anim = None
            self.is_paused = False
            
            # Simulation setup
            self.time_steps = []
            self.lwc_history = []
            self.concentration_history = []
            self.rain_rate_history = []
            
            # Initialize timer
            self.timer = QTimer(self)
            self.timer.timeout.connect(self.start_animation)
            
            # Initialize plot and console
            self.append_log("Cloud Microphysics Window initialized")
            self.update_plot(0)
            
            logging.debug("CloudMicroPhysicsWindow initialization complete")
        except Exception as e:
            logging.error(f"CloudMicroPhysicsWindow initialization failed: {str(e)}")
            self.append_log(f"Error: Failed to initialize window: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize cloud microphysics window: {str(e)}")
            raise
    
    def append_log(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.console.append(f"[{timestamp}] {message}")
        self.console.ensureCursorVisible()
        logging.debug(f"Console log: {message}")
    
    def compute_kohler_activation(self, radii, aerosol_concentration, soluble_fraction):
        logging.debug("Computing Köhler activation")
        try:
            # Köhler theory: critical radius for activation (simplified)
            A = 3.3e-5 / self.air_temp  # Curvature term (μm)
            B = 4.3e-6 * soluble_fraction  # Solute term (μm³)
            critical_radii = np.sqrt(3 * B / A)  # Critical radius (μm)
            activated = radii > critical_radii
            self.droplet_concentrations[activated] = aerosol_concentration * 0.5  # Activate 50% of aerosols
            return self.droplet_concentrations
        except Exception as e:
            logging.error(f"Köhler activation computation failed: {str(e)}")
            self.append_log(f"Error: Köhler activation failed: {str(e)}")
            raise
    
    def compute_condensation_growth(self, radii, supersaturation, air_temp):
        logging.debug("Computing condensation growth")
        try:
            # Clausius-Clapeyron for saturation vapor pressure (Pa)
            es = 611.2 * np.exp((self.Lv / self.Rv) * (1/273.15 - 1/air_temp))
            # Condensation growth: dr/dt = G * S / r, with ventilation effect
            G = 1.0 / (self.rho_water * self.Lv / (es * air_temp) + 1.0)  # Simplified
            ventilation = 1 + 0.1 * self.updraft_vel  # Ventilation factor
            dr_dt = G * supersaturation * ventilation / np.clip(radii, 1e-6, None)
            return np.clip(dr_dt, -1e-3, 1e-3)
        except Exception as e:
            logging.error(f"Condensation growth computation failed: {str(e)}")
            self.append_log(f"Error: Condensation growth failed: {str(e)}")
            raise
    
    def compute_collision_coalescence(self, radii, concentrations):
        logging.debug("Computing collision-coalescence")
        try:
            # Stochastic collision-coalescence: simplified increase in radius
            collision_prob = 0.01 * (radii / self.autoconversion_threshold)**2
            collision_prob = np.clip(collision_prob, 0, 0.1)
            collisions = np.random.random(self.n_droplets) < collision_prob
            radii[collisions] *= 1.2  # Increase radius by 20% for collisions
            concentrations[collisions] *= 0.8  # Reduce concentration
            return radii, concentrations
        except Exception as e:
            logging.error(f"Collision-coalescence computation failed: {str(e)}")
            self.append_log(f"Error: Collision-coalescence failed: {str(e)}")
            raise
    
    def compute_autoconversion(self, radii, concentrations):
        logging.debug("Computing autoconversion")
        try:
            # Kessler-type autoconversion: cloud water to rain
            mask = radii > self.autoconversion_threshold
            auto_rate = 0.001 * np.sum(radii[mask] - self.autoconversion_threshold) / self.autoconversion_threshold
            self.rain_rate = auto_rate * np.sum(self.lwc[mask]) * 3600 / self.rho_water  # mm/hr
            concentrations[mask] *= np.exp(-auto_rate * self.dt)
            return concentrations
        except Exception as e:
            logging.error(f"Autoconversion computation failed: {str(e)}")
            self.append_log(f"Error: Autoconversion failed: {str(e)}")
            raise
    
    def compute_evaporation(self, radii, supersaturation):
        logging.debug("Computing evaporation")
        try:
            # Evaporation when supersaturation is negative
            if supersaturation < 0:
                dr_dt = 0.1 * supersaturation / np.clip(radii, 1e-6, None)
                return np.clip(dr_dt, -1e-3, 0)
            return np.zeros_like(radii)
        except Exception as e:
            logging.error(f"Evaporation computation failed: {str(e)}")
            self.append_log(f"Error: Evaporation failed: {str(e)}")
            raise
    
    def compute_wind_shear_effect(self, radii, concentrations):
        logging.debug("Computing wind shear effect")
        try:
            # Simplified wind shear: perturb droplet positions
            shear_factor = 1 + 0.1 * self.wind_shear * radii / 10.0
            concentrations *= np.clip(shear_factor, 0.5, 1.5)
            return concentrations
        except Exception as e:
            logging.error(f"Wind shear effect computation failed: {str(e)}")
            self.append_log(f"Error: Wind shear effect failed: {str(e)}")
            raise
    
    def update_supersaturation(self):
        logging.debug("Updating supersaturation")
        try:
            # Dynamic supersaturation: cooling from updraft and depletion from condensation
            cooling_rate = -0.01 * self.updraft_vel  # K/s
            self.air_temp += cooling_rate * self.dt
            es = 611.2 * np.exp((self.Lv / self.Rv) * (1/273.15 - 1/self.air_temp))
            e = es * (1 + self.supersaturation)
            # Depletion due to condensation
            condensation_rate = np.sum(self.lwc) * 1e-3 / self.rho_water  # kg/kg
            self.supersaturation -= condensation_rate * 0.1
            self.supersaturation += (self.Lv / (self.c_p * self.air_temp)) * cooling_rate * self.dt
            self.supersaturation = np.clip(self.supersaturation, -0.01, 0.01)
        except Exception as e:
            logging.error(f"Supersaturation update failed: {str(e)}")
            self.append_log(f"Error: Supersaturation update failed: {str(e)}")
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
        self.append_log(f"Animation speed set to {value} ms")
    
    def start_simulation(self):
        logging.debug("Starting cloud microphysics simulation")
        try:
            self.updraft_vel = float(self.inputs["Updraft Velocity (m/s)"].text())
            self.supersaturation = float(self.inputs["Initial Supersaturation (%)"].text()) / 100
            self.initial_radius = float(self.inputs["Initial Droplet Radius (μm)"].text())
            self.initial_concentration = float(self.inputs["Initial Droplet Concentration (cm⁻³)"].text())
            self.aerosol_concentration = float(self.inputs["Aerosol Concentration (cm⁻³)"].text())
            self.aerosol_soluble_fraction = float(self.inputs["Aerosol Soluble Fraction"].text())
            self.air_temp = float(self.inputs["Air Temperature (K)"].text())
            self.wind_speed = float(self.inputs["Wind Speed (m/s)"].text())
            self.wind_angle = float(self.inputs["Wind Angle (degrees)"].text())
            self.autoconversion_threshold = float(self.inputs["Autoconversion Threshold (μm)"].text())
            self.wind_shear = float(self.inputs["Wind Shear (s⁻¹)"].text())
            
            if (self.updraft_vel < 0 or self.supersaturation < 0 or self.supersaturation > 0.01 or 
                self.initial_radius <= 0 or self.initial_concentration <= 0 or 
                self.aerosol_concentration <= 0 or self.aerosol_soluble_fraction < 0 or 
                self.aerosol_soluble_fraction > 1 or self.air_temp < 200 or self.air_temp > 350 or 
                self.wind_speed < 0 or self.autoconversion_threshold <= 0 or self.wind_shear < 0):
                raise ValueError("Invalid parameters: ensure non-negative values, supersaturation ≤ 1%, soluble fraction 0–1, temperature 200–350 K.")
            
            # Reset simulation data
            self.current_step = 0
            self.time_steps = []
            self.lwc_history = []
            self.concentration_history = []
            self.rain_rate_history = []
            self.droplet_radii = np.random.normal(self.initial_radius, self.initial_radius * 0.1, self.n_droplets)
            self.droplet_concentrations = np.full(self.n_droplets, self.initial_concentration)
            self.lwc = (4/3) * np.pi * (self.droplet_radii * 1e-6)**3 * 1e6 * self.droplet_concentrations
            self.z = np.zeros(self.n_droplets)  # Reset height
            self.rain_rate = 0.0
            
            self.console.clear()
            self.append_log("Simulation started")
            self.append_log(f"Initial parameters: Updraft={self.updraft_vel:.2f} m/s, "
                          f"Supersaturation={self.supersaturation*100:.2f}%, "
                          f"Radius={self.initial_radius:.2f} μm, "
                          f"Concentration={self.initial_concentration:.2f} cm⁻³, "
                          f"Aerosol={self.aerosol_concentration:.2f} cm⁻³")
            
            self.timer.start(self.speed_slider.value())
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(True)
            self.resume_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            self.is_paused = False
            self.status_label.setText("Status: Running")
            logging.debug("Cloud microphysics simulation started")
        except ValueError as e:
            logging.error(f"Start cloud microphysics failed: {str(e)}")
            self.append_log(f"Error: Invalid input: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            logging.error(f"Start cloud microphysics failed: {str(e)}")
            self.append_log(f"Error: Failed to start simulation: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to start cloud microphysics: {str(e)}")
    
    def pause_simulation(self):
        logging.debug("Pausing cloud microphysics simulation")
        try:
            if self.anim is not None:
                self.anim.event_source.stop()
                self.is_paused = True
            self.timer.stop()
            self.pause_button.setEnabled(False)
            self.resume_button.setEnabled(True)
            self.stop_button.setEnabled(True)
            self.status_label.setText("Status: Paused")
            self.append_log("Simulation paused")
            logging.debug("Cloud microphysics simulation paused")
        except Exception as e:
            logging.error(f"Pause cloud microphysics failed: {str(e)}")
            self.append_log(f"Error: Failed to pause simulation: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to pause cloud microphysics: {str(e)}")
    
    def resume_simulation(self):
        logging.debug("Resuming cloud microphysics simulation")
        try:
            if self.is_paused:
                self.anim = FuncAnimation(self.figure, self.update_plot, interval=self.speed_slider.value(), 
                                       blit=False, cache_frame_data=False)
                self.timer.start(self.speed_slider.value())
                self.is_paused = False
            self.pause_button.setEnabled(True)
            self.resume_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            self.status_label.setText("Status: Running")
            self.append_log("Simulation resumed")
            logging.debug("Cloud microphysics simulation resumed")
        except Exception as e:
            logging.error(f"Resume cloud microphysics failed: {str(e)}")
            self.append_log(f"Error: Failed to resume simulation: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to resume cloud microphysics: {str(e)}")
    
    def stop_simulation(self):
        logging.debug("Stopping cloud microphysics simulation")
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
            self.status_label.setText("Status: Stopped")
            self.append_log("Simulation stopped")
            logging.debug("Cloud microphysics simulation stopped")
        except Exception as e:
            logging.error(f"Stop cloud microphysics failed: {str(e)}")
            self.append_log(f"Error: Failed to stop simulation: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to stop cloud microphysics: {str(e)}")
    
    def start_animation(self):
        if self.anim is None and not self.is_paused:
            self.anim = FuncAnimation(self.figure, self.update_plot, interval=self.speed_slider.value(), 
                                   blit=False, cache_frame_data=False)
            self.canvas.draw()
    
    def update_plot(self, frame):
        logging.debug(f"Updating cloud microphysics plot at step {self.current_step}")
        try:
            if self.current_step >= int(self.total_time / self.dt):
                self.stop_simulation()
                return
            
            # Update droplet properties
            self.droplet_concentrations = self.compute_kohler_activation(self.droplet_radii, 
                                                                      self.aerosol_concentration, 
                                                                      self.aerosol_soluble_fraction)
            dr_dt = self.compute_condensation_growth(self.droplet_radii, self.supersaturation, self.air_temp)
            dr_dt += self.compute_evaporation(self.droplet_radii, self.supersaturation)
            self.droplet_radii += dr_dt * self.dt
            self.droplet_radii, self.droplet_concentrations = self.compute_collision_coalescence(
                self.droplet_radii, self.droplet_concentrations)
            self.droplet_concentrations = self.compute_autoconversion(self.droplet_radii, self.droplet_concentrations)
            self.droplet_concentrations = self.compute_wind_shear_effect(self.droplet_radii, self.droplet_concentrations)
            self.lwc = (4/3) * np.pi * (self.droplet_radii * 1e-6)**3 * 1e6 * self.droplet_concentrations
            self.z += self.updraft_vel * self.dt  # Update height
            self.update_supersaturation()
            
            # Update status and console
            mean_radius = np.mean(self.droplet_radii)
            mean_lwc = np.mean(self.lwc)
            mean_concentration = np.mean(self.droplet_concentrations)
            mean_height = np.mean(self.z)
            self.status_label.setText(f"Status: Running (Mean Radius: {mean_radius:.2f} μm, "
                                    f"LWC: {mean_lwc:.2f} g/m³, Rain Rate: {self.rain_rate:.2f} mm/hr)")
            self.append_log(f"Step {self.current_step}: Time={self.current_step * self.dt:.2f} s, "
                          f"Mean Radius={mean_radius:.2f} μm, LWC={mean_lwc:.2f} g/m³, "
                          f"Concentration={mean_concentration:.2f} cm⁻³, Rain Rate={self.rain_rate:.2f} mm/hr, "
                          f"Supersaturation={self.supersaturation*100:.2f}%, Height={mean_height:.2f} m")
            
            # Wind stress vector
            angle = self.wind_angle + 0.1 * self.current_step * self.dt
            u_air = self.wind_speed * np.cos(np.radians(angle))
            v_air = self.wind_speed * np.sin(np.radians(angle))
            tau_x = self.rho_air * self.drag_coeff * self.wind_speed * u_air
            tau_y = self.rho_air * self.drag_coeff * self.wind_speed * v_air
            
            # Store history for time series
            self.time_steps.append(self.current_step * self.dt)
            self.lwc_history.append(mean_lwc)
            self.concentration_history.append(mean_concentration)
            self.rain_rate_history.append(self.rain_rate)
            
            # Update plots
            self.ax1.clear()
            self.ax2.clear()
            self.ax3.clear()
            
            plot_style = self.plot_selector.currentText()
            
            # Scatterplot or histogram
            if plot_style == "Color by LWC":
                scatter = self.ax1.scatter(self.droplet_radii, self.droplet_concentrations, 
                                       c=self.lwc, cmap='Blues', s=self.droplet_radii*2, alpha=0.6)
                self.figure.colorbar(scatter, ax=self.ax1, label="LWC (g/m³)")
                self.ax1.set_xlabel("Droplet Radius (μm)")
                self.ax1.set_ylabel("Droplet Concentration (cm⁻³)")
            elif plot_style == "Color by Droplet Size":
                scatter = self.ax1.scatter(self.droplet_radii, self.droplet_concentrations, 
                                       c=self.droplet_radii, cmap='viridis', s=self.droplet_radii*2, alpha=0.6)
                self.figure.colorbar(scatter, ax=self.ax1, label="Droplet Radius (μm)")
                self.ax1.set_xlabel("Droplet Radius (μm)")
                self.ax1.set_ylabel("Droplet Concentration (cm⁻³)")
            elif plot_style == "Droplet Size Histogram":
                self.ax1.hist(self.droplet_radii, bins=20, range=(0, 100), color='blue', alpha=0.6)
                self.ax1.set_xlabel("Droplet Radius (μm)")
                self.ax1.set_ylabel("Frequency")
            else:  # Vertical Profile
                scatter = self.ax1.scatter(self.droplet_radii, self.z, 
                                       c=self.lwc, cmap='Blues', s=self.droplet_radii*2, alpha=0.6)
                self.figure.colorbar(scatter, ax=self.ax1, label="LWC (g/m³)")
                self.ax1.set_xlabel("Droplet Radius (μm)")
                self.ax1.set_ylabel("Height (m)")
            
            self.ax1.set_title(f"Cloud Microphysics at t={self.current_step * self.dt:.2f}s")
            self.ax1.set_xlim(0, 100)
            self.ax1.set_ylim(0, 1000 if plot_style != "Vertical Profile" else max(100, np.max(self.z)))
            
            # Wind stress vector field (not for histogram)
            if plot_style != "Droplet Size Histogram":
                x = np.linspace(0, 100, 5)
                y = np.linspace(0, 1000 if plot_style != "Vertical Profile" else np.max(self.z), 5)
                X, Y = np.meshgrid(x, y)
                U = np.full_like(X, tau_x)
                V = np.full_like(Y, tau_y)
                self.ax1.quiver(X, Y, U, V, color='black', scale=0.5, width=0.002)
            
            # Histogram or vertical profile
            if plot_style != "Droplet Size Histogram":
                self.ax2.hist(self.droplet_radii, bins=20, range=(0, 100), color='blue', alpha=0.6)
                self.ax2.set_xlabel("Droplet Radius (μm)")
                self.ax2.set_ylabel("Frequency")
                self.ax2.set_title("Droplet Size Distribution")
            else:
                self.ax2.scatter(self.droplet_radii, self.z, 
                               c=self.lwc, cmap='Blues', s=self.droplet_radii*2, alpha=0.6)
                self.figure.colorbar(scatter, ax=self.ax2, label="LWC (g/m³)")
                self.ax2.set_xlabel("Droplet Radius (μm)")
                self.ax2.set_ylabel("Height (m)")
                self.ax2.set_title("Vertical Profile")
            
            # Time series plot
            self.ax3.plot(self.time_steps, self.lwc_history, label="Mean LWC (g/m³)")
            self.ax3.plot(self.time_steps, self.concentration_history, label="Mean Concentration (cm⁻³)")
            self.ax3.plot(self.time_steps, self.rain_rate_history, label="Rain Rate (mm/hr)")
            self.ax3.set_xlabel("Time (s)")
            self.ax3.set_ylabel("Values")
            self.ax3.set_title("Time Series")
            self.ax3.legend()
            self.ax3.grid(True)
            
            self.figure.subplots_adjust(left=0.1, right=0.95, wspace=0.3, top=0.95, bottom=0.1, hspace=0.3)
            
            self.current_step += 1
            logging.debug(f"Cloud microphysics plot updated successfully at step {self.current_step}")
            
            return [scatter] if plot_style != "Droplet Size Histogram" else []
        except ValueError as e:
            self.stop_simulation()
            logging.error(f"Cloud microphysics plot update failed: {str(e)}")
            self.append_log(f"Error: Plot update failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            self.stop_simulation()
            logging.error(f"Cloud microphysics plot update failed: {str(e)}")
            self.append_log(f"Error: Plot update failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to update cloud microphysics plot: {str(e)}")
    
    def closeEvent(self, event):
        logging.debug("Closing cloud microphysics window")
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
            self.status_label.setText("Status: Closed")
            self.append_log("Cloud Microphysics Window closed")
            event.accept()
        except Exception as e:
            logging.error(f"Close cloud microphysics window failed: {str(e)}")
            self.append_log(f"Error: Failed to close window: {str(e)}")
            event.accept()