import numpy as np
import logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QMessageBox, QGroupBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
import matplotlib

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class SurfaceLayerPhysicsWindow(QMainWindow):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        logging.debug("Initializing SurfaceLayerPhysicsWindow")
        try:
            self.model = model
            self.setWindowTitle("Surface Layer Physics Analysis")
            self.setGeometry(150, 150, 900, 600)
            
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
            """)
            
            # Main widget and layout
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_layout = QHBoxLayout()
            main_widget.setLayout(main_layout)
            
            # Left side: Control panel
            left_layout = QVBoxLayout()
            control_group = QGroupBox("Surface Layer Parameters")
            form_layout = QFormLayout()
            form_layout.setHorizontalSpacing(10)
            
            self.inputs = {}
            params = [
                ("Drag Coefficient", str(self.model.coupling.drag_coeff)),
                ("Wind Speed (m/s)", str(self.model.coupling.wind_speed)),
                ("Sensible Heat Coeff (W/m²/K)", "10.0"),
                ("Latent Heat Coeff (W/m²)", "20.0")
            ]
            
            for label, default in params:
                edit = QLineEdit(default)
                edit.setFont(QFont("Consolas", 10))
                form_layout.addRow(QLabel(label), edit)
                self.inputs[label] = edit
            
            # Control buttons for dynamic simulation
            button_layout = QHBoxLayout()
            self.start_button = QPushButton("Start Simulation")
            self.start_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.start_button.clicked.connect(self.start_simulation)
            
            self.stop_button = QPushButton("Stop Simulation")
            self.stop_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.stop_button.clicked.connect(self.stop_simulation)
            self.stop_button.setEnabled(False)
            
            button_layout.addWidget(self.start_button)
            button_layout.addWidget(self.stop_button)
            
            form_layout.addRow(button_layout)
            control_group.setLayout(form_layout)
            left_layout.addWidget(control_group)
            left_layout.addStretch()
            
            # Right side: Plot widget
            self.figure = Figure(facecolor="#C0C0C0")
            self.canvas = FigureCanvas(self.figure)
            matplotlib.rc("axes", edgecolor="black", labelcolor="black", facecolor="#FFFFFF")
            matplotlib.rc("xtick", color="black")
            matplotlib.rc("ytick", color="black")
            matplotlib.rc("text", color="black")
            matplotlib.rc("figure", facecolor="#C0C0C0")
            matplotlib.rc("font", family="Consolas", size=10)
            
            # Add layouts to main layout
            main_layout.addLayout(left_layout, 1)
            main_layout.addWidget(self.canvas, 3)
            
            # Initialize plot elements
            self.ax1 = self.figure.add_subplot(211)
            self.ax2 = self.figure.add_subplot(212)
            self.current_step = 0
            self.anim = None
            
            # Simulation setup
            self.time_steps = []
            self.heat_fluxes = []
            self.wind_stresses = []
            self.dt = 1800.0  # Match main simulation time step
            self.ny, self.nx = self.model.ocean_temps.shape
            self.surface_currents = np.zeros((self.ny, self.nx))
            
            # Initialize plot
            self.update_plot(0)
            
            # Timer for animation
            self.timer = QTimer()
            self.timer.timeout.connect(self.start_animation)
            
            logging.debug("SurfaceLayerPhysicsWindow initialization complete")
        except Exception as e:
            logging.error(f"SurfaceLayerPhysicsWindow initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize surface layer physics window: {str(e)}")
            raise
    
    def start_simulation(self):
        logging.debug("Starting surface layer physics simulation")
        try:
            drag_coeff = float(self.inputs["Drag Coefficient"].text())
            wind_speed = float(self.inputs["Wind Speed (m/s)"].text())
            sensible_heat_coeff = float(self.inputs["Sensible Heat Coeff (W/m²/K)"].text())
            latent_heat_coeff = float(self.inputs["Latent Heat Coeff (W/m²)"].text())
            
            if drag_coeff <= 0 or wind_speed < 0 or sensible_heat_coeff < 0 or latent_heat_coeff < 0:
                raise ValueError("Parameters must be non-negative; drag coefficient must be positive.")
            
            self.model.coupling.drag_coeff = drag_coeff
            self.model.coupling.wind_speed = wind_speed
            self.current_step = 0
            self.time_steps = []
            self.heat_fluxes = []
            self.wind_stresses = []
            self.surface_currents = np.zeros((self.ny, self.nx))
            self.timer.start(100)
            self.start_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            logging.debug("Surface layer physics simulation started")
        except ValueError as e:
            logging.error(f"Start surface layer physics failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            logging.error(f"Start surface layer physics failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to start surface layer physics: {str(e)}")
    
    def stop_simulation(self):
        logging.debug("Stopping surface layer physics simulation")
        try:
            self.timer.stop()
            if self.anim is not None:
                self.anim.event_source.stop()
                self.anim = None
            self.start_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            logging.debug("Surface layer physics simulation stopped")
        except Exception as e:
            logging.error(f"Stop surface layer physics failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to stop surface layer physics: {str(e)}")
    
    def start_animation(self):
        if self.anim is None:
            self.anim = FuncAnimation(self.figure, self.update_plot, interval=100, 
                                   blit=False, cache_frame_data=False)
            self.canvas.draw()
    
    def update_plot(self, frame):
        logging.debug(f"Updating surface layer physics plot at step {self.current_step}")
        try:
            if not self.model or self.current_step >= int(self.model.total_time / self.model.dt):
                self.stop_simulation()
                return
            
            drag_coeff = float(self.inputs["Drag Coefficient"].text())
            wind_speed = float(self.inputs["Wind Speed (m/s)"].text())
            sensible_heat_coeff = float(self.inputs["Sensible Heat Coeff (W/m²/K)"].text())
            latent_heat_coeff = float(self.inputs["Latent Heat Coeff (W/m²)"].text())
            
            if drag_coeff <= 0 or wind_speed < 0 or sensible_heat_coeff < 0 or latent_heat_coeff < 0:
                self.stop_simulation()
                raise ValueError("Parameters must be non-negative; drag coefficient must be positive.")
            
            # Compute surface layer physics
            ocean_temps = self.model.ocean_temps
            atm_temps = self.model.atm_temps
            rho_air = 1.225  # Air density (kg/m³)
            Cp_air = 1005.0  # Specific heat of air (J/kg/K)
            L_vap = 2.5e6    # Latent heat of vaporization (J/kg)
            
            # Sensible heat flux: Q_s = rho_air * Cp_air * C_h * U * (T_ocean - T_atm)
            sensible_heat = rho_air * Cp_air * sensible_heat_coeff * wind_speed * (ocean_temps - atm_temps)
            
            # Latent heat flux: Q_l = rho_air * L_vap * C_e * U * q
            # Simplified: assume q (specific humidity difference) is constant for demo
            q = 0.001  # Placeholder for humidity difference
            latent_heat = rho_air * L_vap * latent_heat_coeff * wind_speed * q
            
            # Total heat flux
            total_heat_flux = sensible_heat + latent_heat
            
            # Wind stress: tau = rho_air * C_d * U^2
            wind_stress = rho_air * drag_coeff * wind_speed**2
            
            # Update surface currents (simplified: proportional to wind stress)
            mixing_coeff = self.model.coupling.mixing_coeff
            self.surface_currents += self.dt * (wind_stress / 1000.0) * np.random.normal(0, mixing_coeff, (self.ny, self.nx))
            self.surface_currents = np.clip(self.surface_currents, -0.5, 0.5)  # Prevent unrealistic values
            
            # Store for time series
            self.time_steps.append(self.current_step * self.dt)
            self.heat_fluxes.append(np.mean(total_heat_flux))
            self.wind_stresses.append(wind_stress)
            
            # Update plots
            self.ax1.clear()
            self.ax2.clear()
            
            # Heatmap of surface currents
            im = self.ax1.imshow(self.surface_currents, cmap="coolwarm", origin="lower", vmin=-0.5, vmax=0.5)
            self.ax1.set_title(f"Surface Currents (m/s) at t={self.current_step * self.dt:.2f}s")
            self.ax1.set_xlabel("X")
            self.ax1.set_ylabel("Y")
            if not hasattr(self.ax1, 'colorbar'):
                self.figure.colorbar(im, ax=self.ax1, label="Current Speed (m/s)")
            
            # Time series of heat flux and wind stress
            self.ax2.plot(self.time_steps, self.heat_fluxes, label="Mean Heat Flux (W/m²)", color="blue")
            self.ax2.plot(self.time_steps, self.wind_stresses, label="Wind Stress (N/m²)", color="red")
            self.ax2.set_title("Surface Layer Dynamics Over Time")
            self.ax2.set_xlabel("Time (s)")
            self.ax2.set_ylabel("Flux/Stress")
            self.ax2.legend()
            self.ax2.grid(True, color="gray", linestyle="--", alpha=0.5)
            
            self.figure.subplots_adjust(left=0.1, right=0.95, hspace=0.3, top=0.9, bottom=0.1)
            
            self.current_step += 1
            logging.debug(f"Surface layer physics plot updated at step {self.current_step}")
            
            return [im]
        except ValueError as e:
            self.stop_simulation()
            logging.error(f"Surface layer physics plot update failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            self.stop_simulation()
            logging.error(f"Surface layer physics plot update failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to update surface layer physics plot: {str(e)}")
    
    def closeEvent(self, event):
        logging.debug("Closing surface layer physics window")
        self.timer.stop()
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
        self.start_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        event.accept()