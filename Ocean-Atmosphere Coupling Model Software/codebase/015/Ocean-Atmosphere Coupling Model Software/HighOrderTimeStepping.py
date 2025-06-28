import numpy as np
import logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QMessageBox, QGroupBox, QComboBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
import matplotlib

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class HighOrderTimeSteppingWindow(QMainWindow):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        logging.debug("Initializing HighOrderTimeSteppingWindow")
        try:
            self.model = model
            self.setWindowTitle("High-Order Time Stepping Analysis")
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
            """)
            
            # Main widget and layout
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_layout = QHBoxLayout()
            main_widget.setLayout(main_layout)
            
            # Left side: Control panel
            left_layout = QVBoxLayout()
            control_group = QGroupBox("Time Stepping Parameters")
            form_layout = QFormLayout()
            form_layout.setHorizontalSpacing(10)
            
            self.inputs = {}
            params = [
                ("Time Step (s)", "1800.0"),
                ("Drag Coefficient", str(self.model.coupling.drag_coeff)),
                ("Sensible Heat Coeff (W/m²/K)", "10.0"),
                ("Boundary Layer Depth (m)", "50.0"),
                ("KPP Mixing Coeff (m²/s)", "0.01")
            ]
            
            for label, default in params:
                edit = QLineEdit(default)
                edit.setFont(QFont("Consolas", 10))
                form_layout.addRow(QLabel(label), edit)
                self.inputs[label] = edit
            
            # Scheme selection
            self.scheme_selector = QComboBox()
            self.scheme_selector.addItems(["Bulk Scheme", "KPP Scheme"])
            self.scheme_selector.setFont(QFont("Consolas", 10))
            form_layout.addRow(QLabel("Scheme"), self.scheme_selector)
            
            # Time stepping method selection
            self.method_selector = QComboBox()
            self.method_selector.addItems(["Euler", "RK4"])
            self.method_selector.setFont(QFont("Consolas", 10))
            form_layout.addRow(QLabel("Time Stepping Method"), self.method_selector)
            
            # Control buttons
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
            matplotlib.rc("font", family="Consolas", size=10)
            
            # Add layouts to main layout
            main_layout.addLayout(left_layout, 1)
            main_layout.addWidget(self.canvas, 3)
            
            # Initialize plot elements
            self.ax1 = self.figure.add_subplot(121)
            self.ax2 = self.figure.add_subplot(122)
            self.current_step = 0
            self.anim = None
            
            # Simulation setup
            self.time_steps = []
            self.euler_values = []
            self.rk4_values = []
            self.dt = 1800.0
            self.ny, self.nx = self.model.ocean_temps.shape
            self.bulk_heat_flux_euler = np.zeros((self.ny, self.nx))
            self.kpp_diffusivity_euler = np.zeros((self.ny, self.nx))
            self.bulk_heat_flux_rk4 = np.zeros((self.ny, self.nx))
            self.kpp_diffusivity_rk4 = np.zeros((self.ny, self.nx))
            
            # Initialize plot
            self.update_plot(0)
            
            # Timer for animation
            self.timer = QTimer()
            self.timer.timeout.connect(self.start_animation)
            
            logging.debug("HighOrderTimeSteppingWindow initialization complete")
        except Exception as e:
            logging.error(f"HighOrderTimeSteppingWindow initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize high-order time stepping window: {str(e)}")
            raise
    
    def start_simulation(self):
        logging.debug("Starting high-order time stepping simulation")
        try:
            dt = float(self.inputs["Time Step (s)"].text())
            drag_coeff = float(self.inputs["Drag Coefficient"].text())
            sensible_heat_coeff = float(self.inputs["Sensible Heat Coeff (W/m²/K)"].text())
            boundary_layer_depth = float(self.inputs["Boundary Layer Depth (m)"].text())
            kpp_mixing_coeff = float(self.inputs["KPP Mixing Coeff (m²/s)"].text())
            
            if dt <= 0 or drag_coeff <= 0 or sensible_heat_coeff < 0 or boundary_layer_depth <= 0 or kpp_mixing_coeff <= 0:
                raise ValueError("Parameters must be positive; time step, drag coefficient, and mixing coefficient must be non-zero.")
            
            self.dt = dt
            self.model.coupling.drag_coeff = drag_coeff
            self.current_step = 0
            self.time_steps = []
            self.euler_values = []
            self.rk4_values = []
            self.bulk_heat_flux_euler = np.zeros((self.ny, self.nx))
            self.kpp_diffusivity_euler = np.zeros((self.ny, self.nx))
            self.bulk_heat_flux_rk4 = np.zeros((self.ny, self.nx))
            self.kpp_diffusivity_rk4 = np.zeros((self.ny, self.nx))
            self.timer.start(100)
            self.start_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            logging.debug("High-order time stepping simulation started")
        except ValueError as e:
            logging.error(f"Start high-order time stepping failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            logging.error(f"Start high-order time stepping failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to start high-order time stepping: {str(e)}")
    
    def stop_simulation(self):
        logging.debug("Stopping high-order time stepping simulation")
        try:
            self.timer.stop()
            if self.anim is not None:
                self.anim.event_source.stop()
                self.anim = None
            self.start_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            logging.debug("High-order time stepping simulation stopped")
        except Exception as e:
            logging.error(f"Stop high-order time stepping failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to stop high-order time stepping: {str(e)}")
    
    def start_animation(self):
        if self.anim is None:
            self.anim = FuncAnimation(self.figure, self.update_plot, interval=100, 
                                   blit=False, cache_frame_data=False)
            self.canvas.draw()
    
    def update_plot(self, frame):
        logging.debug(f"Updating high-order time stepping plot at step {self.current_step}")
        try:
            if not self.model or self.current_step >= int(self.model.total_time / self.dt):
                self.stop_simulation()
                return
            
            dt = float(self.inputs["Time Step (s)"].text())
            sensible_heat_coeff = float(self.inputs["Sensible Heat Coeff (W/m²/K)"].text())
            boundary_layer_depth = float(self.inputs["Boundary Layer Depth (m)"].text())
            kpp_mixing_coeff = float(self.inputs["KPP Mixing Coeff (m²/s)"].text())
            scheme = self.scheme_selector.currentText()
            method = self.method_selector.currentText()
            
            if dt <= 0 or sensible_heat_coeff < 0 or boundary_layer_depth <= 0 or kpp_mixing_coeff <= 0:
                self.stop_simulation()
                raise ValueError("Parameters must be positive; time step and mixing coefficient must be non-zero.")
            
            ocean_temps = self.model.ocean_temps
            atm_temps = self.model.atm_temps
            wind_speed = self.model.coupling.wind_speed
            rho_air = 1.225  # Air density (kg/m³)
            Cp_air = 1005.0  # Specific heat of air (J/kg/K)
            
            # Euler method
            bulk_heat_flux_euler = rho_air * Cp_air * sensible_heat_coeff * wind_speed * (ocean_temps - atm_temps)
            z = np.linspace(0, boundary_layer_depth, self.ny)
            K_euler = kpp_mixing_coeff * (1 - z[:, np.newaxis] / boundary_layer_depth)**2
            K_euler = np.clip(K_euler, 0, kpp_mixing_coeff)
            temp_gradient = np.gradient(ocean_temps, axis=0) / self.model.dy
            kpp_heat_flux_euler = -K_euler * temp_gradient
            
            # RK4 method
            def compute_flux_derivative(ocean_temps, atm_temps, is_kpp=False):
                if is_kpp:
                    K = kpp_mixing_coeff * (1 - z[:, np.newaxis] / boundary_layer_depth)**2
                    K = np.clip(K, 0, kpp_mixing_coeff)
                    temp_gradient = np.gradient(ocean_temps, axis=0) / self.model.dy
                    return -K * temp_gradient
                else:
                    return rho_air * Cp_air * sensible_heat_coeff * wind_speed * (ocean_temps - atm_temps)
            
            # RK4 integration
            k1 = compute_flux_derivative(ocean_temps, atm_temps, is_kpp=(scheme == "KPP Scheme"))
            k2 = compute_flux_derivative(ocean_temps + 0.5 * dt * k1 / sensible_heat_coeff if not (scheme == "KPP Scheme") else ocean_temps,
                                       atm_temps + 0.5 * dt * k1 / sensible_heat_coeff if not (scheme == "KPP Scheme") else atm_temps,
                                       is_kpp=(scheme == "KPP Scheme"))
            k3 = compute_flux_derivative(ocean_temps + 0.5 * dt * k2 / sensible_heat_coeff if not (scheme == "KPP Scheme") else ocean_temps,
                                       atm_temps + 0.5 * dt * k2 / sensible_heat_coeff if not (scheme == "KPP Scheme") else atm_temps,
                                       is_kpp=(scheme == "KPP Scheme"))
            k4 = compute_flux_derivative(ocean_temps + dt * k3 / sensible_heat_coeff if not (scheme == "KPP Scheme") else ocean_temps,
                                       atm_temps + dt * k3 / sensible_heat_coeff if not (scheme == "KPP Scheme") else atm_temps,
                                       is_kpp=(scheme == "KPP Scheme"))
            if scheme == "Bulk Scheme":
                bulk_heat_flux_rk4 = bulk_heat_flux_euler + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
            else:
                kpp_heat_flux_rk4 = kpp_heat_flux_euler + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
            
            # Update fields
            if scheme == "Bulk Scheme":
                self.bulk_heat_flux_euler = bulk_heat_flux_euler
                self.bulk_heat_flux_rk4 = bulk_heat_flux_rk4
            else:
                self.kpp_diffusivity_euler = K_euler
                self.kpp_diffusivity_rk4 = K_euler  # Diffusivity is static, but heat flux evolves
            
            # Store mean values for line plot
            self.time_steps.append(self.current_step * self.dt)
            if scheme == "Bulk Scheme":
                self.euler_values.append(np.mean(self.bulk_heat_flux_euler))
                self.rk4_values.append(np.mean(self.bulk_heat_flux_rk4))
            else:
                self.euler_values.append(np.mean(kpp_heat_flux_euler))
                self.rk4_values.append(np.mean(kpp_heat_flux_rk4))
            
            # Update plots
            self.ax1.clear()
            self.ax2.clear()
            
            # Heatmap with annotations
            x = np.arange(self.nx)
            y = np.arange(self.ny)
            if scheme == "Bulk Scheme":
                data = self.bulk_heat_flux_euler if method == "Euler" else self.bulk_heat_flux_rk4
                title = f"Bulk Heat Flux ({method}, W/m²)\nat t={self.current_step * self.dt:.2f}s"
                cmap = "coolwarm"
                label = "Heat Flux (W/m²)"
            else:
                data = self.kpp_diffusivity_euler if method == "Euler" else self.kpp_diffusivity_rk4
                title = f"KPP Diffusivity ({method}, m²/s)\nat t={self.current_step * self.dt:.2f}s"
                cmap = "viridis"
                label = "Diffusivity (m²/s)"
            
            im = self.ax1.imshow(data, cmap=cmap, origin="lower")
            self.ax1.set_title(title)
            self.ax1.set_xlabel("X")
            self.ax1.set_ylabel("Y")
            if not hasattr(self.ax1, 'colorbar'):
                self.figure.colorbar(im, ax=self.ax1, label=label)
            
            # Add annotations
            for i in range(0, self.ny, max(1, self.ny // 5)):
                for j in range(0, self.nx, max(1, self.nx // 5)):
                    self.ax1.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center", color="black" if data[i, j] < np.mean(data) else "white")
            
            # Line plot: Mean values over time
            self.ax2.plot(self.time_steps, self.euler_values, label="Euler", color="blue")
            self.ax2.plot(self.time_steps, self.rk4_values, label="RK4", color="red")
            self.ax2.set_title(f"Mean {label} Over Time")
            self.ax2.set_xlabel("Time (s)")
            self.ax2.set_ylabel(label)
            self.ax2.legend()
            self.ax2.grid(True, color="gray", linestyle="--", alpha=0.5)
            
            self.figure.subplots_adjust(left=0.1, right=0.95, wspace=0.3, top=0.9, bottom=0.1)
            
            self.current_step += 1
            logging.debug(f"High-order time stepping plot updated at step {self.current_step}")
            
            return [im]
        except ValueError as e:
            self.stop_simulation()
            logging.error(f"High-order time stepping plot update failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            self.stop_simulation()
            logging.error(f"High-order time stepping plot update failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to update high-order time stepping plot: {str(e)}")
    
    def closeEvent(self, event):
        logging.debug("Closing high-order time stepping window")
        self.timer.stop()
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
        self.start_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        event.accept()