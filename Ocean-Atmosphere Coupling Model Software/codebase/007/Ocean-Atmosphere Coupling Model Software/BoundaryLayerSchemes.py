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

class BoundaryLayerSchemesWindow(QMainWindow):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        logging.debug("Initializing BoundaryLayerSchemesWindow")
        try:
            self.model = model
            self.setWindowTitle("Boundary Layer Schemes Analysis")
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
            control_group = QGroupBox("Boundary Layer Parameters")
            form_layout = QFormLayout()
            form_layout.setHorizontalSpacing(10)
            
            self.inputs = {}
            params = [
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
            self.bulk_fluxes = []
            self.kpp_fluxes = []
            self.dt = 1800.0  # Match main simulation time step
            self.ny, self.nx = self.model.ocean_temps.shape
            self.bulk_heat_flux = np.zeros((self.ny, self.nx))
            self.kpp_diffusivity = np.zeros((self.ny, self.nx))
            
            # Initialize plot
            self.update_plot(0)
            
            # Timer for animation
            self.timer = QTimer()
            self.timer.timeout.connect(self.start_animation)
            
            logging.debug("BoundaryLayerSchemesWindow initialization complete")
        except Exception as e:
            logging.error(f"BoundaryLayerSchemesWindow initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize boundary layer schemes window: {str(e)}")
            raise
    
    def start_simulation(self):
        logging.debug("Starting boundary layer schemes simulation")
        try:
            drag_coeff = float(self.inputs["Drag Coefficient"].text())
            sensible_heat_coeff = float(self.inputs["Sensible Heat Coeff (W/m²/K)"].text())
            boundary_layer_depth = float(self.inputs["Boundary Layer Depth (m)"].text())
            kpp_mixing_coeff = float(self.inputs["KPP Mixing Coeff (m²/s)"].text())
            
            if drag_coeff <= 0 or sensible_heat_coeff < 0 or boundary_layer_depth <= 0 or kpp_mixing_coeff <= 0:
                raise ValueError("Parameters must be positive; drag coefficient and mixing coefficient must be non-zero.")
            
            self.model.coupling.drag_coeff = drag_coeff
            self.current_step = 0
            self.time_steps = []
            self.bulk_fluxes = []
            self.kpp_fluxes = []
            self.bulk_heat_flux = np.zeros((self.ny, self.nx))
            self.kpp_diffusivity = np.zeros((self.ny, self.nx))
            self.timer.start(100)
            self.start_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            logging.debug("Boundary layer schemes simulation started")
        except ValueError as e:
            logging.error(f"Start boundary layer schemes failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            logging.error(f"Start boundary layer schemes failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to start boundary layer schemes: {str(e)}")
    
    def stop_simulation(self):
        logging.debug("Stopping boundary layer schemes simulation")
        try:
            self.timer.stop()
            if self.anim is not None:
                self.anim.event_source.stop()
                self.anim = None
            self.start_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            logging.debug("Boundary layer schemes simulation stopped")
        except Exception as e:
            logging.error(f"Stop boundary layer schemes failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to stop boundary layer schemes: {str(e)}")
    
    def start_animation(self):
        if self.anim is None:
            self.anim = FuncAnimation(self.figure, self.update_plot, interval=100, 
                                   blit=False, cache_frame_data=False)
            self.canvas.draw()
    
    def update_plot(self, frame):
        logging.debug(f"Updating boundary layer schemes plot at step {self.current_step}")
        try:
            if not self.model or self.current_step >= int(self.model.total_time / self.model.dt):
                self.stop_simulation()
                return
            
            drag_coeff = float(self.inputs["Drag Coefficient"].text())
            sensible_heat_coeff = float(self.inputs["Sensible Heat Coeff (W/m²/K)"].text())
            boundary_layer_depth = float(self.inputs["Boundary Layer Depth (m)"].text())
            kpp_mixing_coeff = float(self.inputs["KPP Mixing Coeff (m²/s)"].text())
            scheme = self.scheme_selector.currentText()
            
            if drag_coeff <= 0 or sensible_heat_coeff < 0 or boundary_layer_depth <= 0 or kpp_mixing_coeff <= 0:
                self.stop_simulation()
                raise ValueError("Parameters must be positive; drag coefficient and mixing coefficient must be non-zero.")
            
            ocean_temps = self.model.ocean_temps
            atm_temps = self.model.atm_temps
            wind_speed = self.model.coupling.wind_speed
            rho_air = 1.225  # Air density (kg/m³)
            Cp_air = 1005.0  # Specific heat of air (J/kg/K)
            
            # Bulk scheme: Sensible heat flux
            bulk_heat_flux = rho_air * Cp_air * sensible_heat_coeff * wind_speed * (ocean_temps - atm_temps)
            
            # KPP scheme: Turbulent diffusivity
            # Simplified KPP: K = k * (1 - z/h)^2 for z < h, where h is boundary layer depth
            z = np.linspace(0, boundary_layer_depth, self.ny)
            K = kpp_mixing_coeff * (1 - z[:, np.newaxis] / boundary_layer_depth)**2
            K = np.clip(K, 0, kpp_mixing_coeff)  # Ensure non-negative
            temp_gradient = np.gradient(ocean_temps, axis=0) / self.model.dy
            kpp_heat_flux = -K * temp_gradient  # Simplified heat flux via KPP
            
            # Update fields based on selected scheme
            if scheme == "Bulk Scheme":
                self.bulk_heat_flux = bulk_heat_flux
            else:  # KPP Scheme
                self.kpp_diffusivity = K
            
            # Store for time series
            self.time_steps.append(self.current_step * self.dt)
            self.bulk_fluxes.append(np.mean(bulk_heat_flux))
            self.kpp_fluxes.append(np.mean(kpp_heat_flux))
            
            # Update plots
            self.ax1.clear()
            self.ax2.clear()
            
            # Contour plot for selected scheme
            x = np.arange(self.nx)
            y = np.arange(self.ny)
            X, Y = np.meshgrid(x, y)
            
            if scheme == "Bulk Scheme":
                im = self.ax1.contourf(X, Y, self.bulk_heat_flux, cmap="coolwarm", levels=20)
                self.ax1.set_title(f"Bulk Scheme Heat Flux (W/m²)\nat t={self.current_step * self.dt:.2f}s")
            else:
                im = self.ax1.contourf(X, Y, self.kpp_diffusivity, cmap="viridis", levels=20)
                self.ax1.set_title(f"KPP Scheme Diffusivity (m²/s)\nat t={self.current_step * self.dt:.2f}s")
            
            self.ax1.set_xlabel("X")
            self.ax1.set_ylabel("Y")
            if not hasattr(self.ax1, 'colorbar'):
                self.figure.colorbar(im, ax=self.ax1, label="Heat Flux (W/m²)" if scheme == "Bulk Scheme" else "Diffusivity (m²/s)")
            
            # Time series of mean fluxes
            self.ax2.plot(self.time_steps, self.bulk_fluxes, label="Bulk Heat Flux (W/m²)", color="blue")
            self.ax2.plot(self.time_steps, self.kpp_fluxes, label="KPP Heat Flux (W/m²)", color="red")
            self.ax2.set_title("Mean Heat Flux Over Time")
            self.ax2.set_xlabel("Time (s)")
            self.ax2.set_ylabel("Heat Flux (W/m²)")
            self.ax2.legend()
            self.ax2.grid(True, color="gray", linestyle="--", alpha=0.5)
            
            self.figure.subplots_adjust(left=0.1, right=0.95, wspace=0.3, top=0.9, bottom=0.1)
            
            self.current_step += 1
            logging.debug(f"Boundary layer schemes plot updated at step {self.current_step}")
            
            return [im]
        except ValueError as e:
            self.stop_simulation()
            logging.error(f"Boundary layer schemes plot update failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            self.stop_simulation()
            logging.error(f"Boundary layer schemes plot update failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to update boundary layer schemes plot: {str(e)}")
    
    def closeEvent(self, event):
        logging.debug("Closing boundary layer schemes window")
        self.timer.stop()
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
        self.start_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        event.accept()