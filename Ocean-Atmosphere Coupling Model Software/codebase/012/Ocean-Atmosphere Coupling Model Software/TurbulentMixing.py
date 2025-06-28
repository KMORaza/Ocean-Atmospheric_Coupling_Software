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

class TurbulentMixingWindow(QMainWindow):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        logging.debug("Initializing TurbulentMixingWindow")
        try:
            self.model = model
            self.setWindowTitle("Turbulent Mixing Analysis")
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
            control_group = QGroupBox("Turbulent Mixing Parameters")
            form_layout = QFormLayout()
            form_layout.setHorizontalSpacing(10)
            
            self.inputs = {}
            params = [
                ("Mixing Coefficient (m²/s)", str(self.model.coupling.mixing_coeff)),
                ("Wind Speed (m/s)", str(self.model.coupling.wind_speed))
            ]
            
            for label, default in params:
                edit = QLineEdit(default)
                edit.setFont(QFont("Consolas", 10))
                form_layout.addRow(QLabel(label), edit)
                self.inputs[label] = edit
            
            # Control buttons for dynamic simulation
            button_layout = QHBoxLayout()
            self.start_button = QPushButton("Start Mixing")
            self.start_button.setFont(QFont("Consolas", 10, QFont.Bold))
            self.start_button.clicked.connect(self.start_simulation)
            
            self.stop_button = QPushButton("Stop Mixing")
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
            
            # Add layouts to main layout
            main_layout.addLayout(left_layout, 1)
            main_layout.addWidget(self.canvas, 3)
            
            # Initialize plot elements
            self.ax1 = self.figure.add_subplot(121)
            self.ax2 = self.figure.add_subplot(122)
            self.current_step = 0
            self.anim = None
            
            # Vertical mixing setup
            self.nz = 50
            self.z = np.linspace(0, -100, self.nz)
            self.dz = self.z[1] - self.z[0]
            self.temp_profile = np.linspace(300, 290, self.nz)
            self.dt = 0.1
            
            # Timer to trigger animation updates
            self.timer = QTimer()
            self.timer.timeout.connect(self.start_animation)
            
            # Initialize plot
            self.update_plot(0)
            
            logging.debug("TurbulentMixingWindow initialization complete")
        except Exception as e:
            logging.error(f"TurbulentMixingWindow initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to initialize turbulent mixing window: {str(e)}")
            raise
    
    def start_simulation(self):
        logging.debug("Starting turbulent mixing simulation")
        try:
            mixing_coeff = float(self.inputs["Mixing Coefficient (m²/s)"].text())
            wind_speed = float(self.inputs["Wind Speed (m/s)"].text())
            
            if mixing_coeff <= 0 or wind_speed < 0:
                raise ValueError("Mixing coefficient must be positive and wind speed must be non-negative.")
            
            self.model.coupling.mixing_coeff = mixing_coeff
            self.model.coupling.wind_speed = wind_speed
            self.current_step = 0
            self.temp_profile = np.linspace(300, 290, self.nz)
            self.timer.start(100)
            self.start_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            logging.debug("Turbulent mixing simulation started")
        except ValueError as e:
            logging.error(f"Start turbulent mixing failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            logging.error(f"Start turbulent mixing failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to start turbulent mixing: {str(e)}")
    
    def stop_simulation(self):
        logging.debug("Stopping turbulent mixing simulation")
        try:
            self.timer.stop()
            if self.anim is not None:
                self.anim.event_source.stop()
                self.anim = None
            self.start_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            logging.debug("Turbulent mixing simulation stopped")
        except Exception as e:
            logging.error(f"Stop turbulent mixing failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to stop turbulent mixing: {str(e)}")
    
    def start_animation(self):
        if self.anim is None:
            self.anim = FuncAnimation(self.figure, self.update_plot, interval=100, 
                                   blit=False, cache_frame_data=False)
            self.canvas.draw()
    
    def update_plot(self, frame):
        logging.debug(f"Updating turbulent mixing plot at step {self.current_step}")
        try:
            if not self.model or self.current_step >= int(self.model.total_time / self.model.dt):
                self.stop_simulation()
                return
            
            mixing_coeff = float(self.inputs["Mixing Coefficient (m²/s)"].text())
            wind_speed = float(self.inputs["Wind Speed (m/s)"].text())
            
            if mixing_coeff <= 0 or wind_speed < 0:
                self.stop_simulation()
                raise ValueError("Mixing coefficient must be positive and wind speed must be non-negative.")
            
            self.model.coupling.mixing_coeff = mixing_coeff
            self.model.coupling.wind_speed = wind_speed
            
            ny, nx = self.model.ocean_temps.shape
            x = np.arange(0, nx)
            y = np.arange(0, ny)
            X, Y = np.meshgrid(x, y)
            x0 = nx / 2 + 10 * np.sin(self.current_step * 0.01 * wind_speed)
            L = 5.0
            T_cold = 290
            T_hot = 300
            front = T_cold + (T_hot - T_cold) * 0.5 * (1 + np.tanh((X - x0) / L))
            ocean_temps = self.model.ocean_temps.copy() * 0.8 + front * 0.2
            
            mixing_field = self.model.coupling.compute_turbulent_mixing(
                ocean_temps, self.model.salinity, self.model.dx, self.model.dy, wind_speed
            )
            
            dy, dx = np.gradient(mixing_field)
            skip = (slice(None, None, 5), slice(None, None, 5))
            X_sub = X[skip]
            Y_sub = Y[skip]
            dx_sub = dx[skip]
            dy_sub = dy[skip]
            
            f = 1e-4
            rho = 1000
            nu = mixing_coeff
            tau = 0.1 * wind_speed**2
            V0 = tau / (rho * np.sqrt(2 * nu * f))
            d = np.sqrt(2 * nu / f)
            u = V0 * np.exp(self.z / d) * np.cos(np.pi / 4 + self.z / d)
            v = V0 * np.exp(self.z / d) * np.sin(np.pi / 4 + self.z / d)
            
            self.temp_profile[0] = ocean_temps[ny//2, nx//2]
            Kv = mixing_coeff * np.exp(self.z / 50)
            d2Tdz2 = (self.temp_profile[2:] - 2 * self.temp_profile[1:-1] + self.temp_profile[:-2]) / self.dz**2
            dTdt = np.zeros(self.nz)
            dTdt[1:-1] = Kv[1:-1] * d2Tdz2
            dTdt[0] = 0
            dTdt[-1] = dTdt[-2]
            self.temp_profile += dTdt * self.dt
            
            self.ax1.clear()
            self.ax2.clear()
            
            im = self.ax1.contourf(X, Y, ocean_temps, cmap="coolwarm", levels=20)
            self.ax1.quiver(X_sub, Y_sub, dx_sub, dy_sub, color="black", scale=20, scale_units="inches")
            self.ax1.set_title(f"Thermal Front\nStep: {self.current_step}")
            self.ax1.set_xlabel("X")
            self.ax1.set_ylabel("Y")
            if not hasattr(self.ax1, 'colorbar'):
                self.figure.colorbar(im, ax=self.ax1, label="Temperature (K)")
            
            self.ax2.plot(u, self.z, label="u (m/s)", color="cyan")
            self.ax2.plot(v, self.z, label="v (m/s)", color="magenta")
            self.ax2.plot(self.temp_profile, self.z, label="Temp (K)", color="yellow")
            self.ax2.set_title("Ekman Spiral & Vertical Mixing")
            self.ax2.set_xlabel("Velocity (m/s) or Temp (K)")
            self.ax2.set_ylabel("Depth (m)")
            self.ax2.legend()
            self.ax2.grid(True, color="gray", linestyle="--", alpha=0.5)
            
            self.figure.subplots_adjust(left=0.1, right=0.95, wspace=0.3, top=0.9, bottom=0.1)
            
            self.current_step += 1
            logging.debug(f"Turbulent mixing plot updated at step {self.current_step}")
            
            return [im]
        except ValueError as e:
            self.stop_simulation()
            logging.error(f"Turbulent mixing plot update failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))
        except Exception as e:
            self.stop_simulation()
            logging.error(f"Turbulent mixing plot update failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to update turbulent mixing plot: {str(e)}")
    
    def closeEvent(self, event):
        logging.debug("Closing turbulent mixing window")
        self.timer.stop()
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
        self.start_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        event.accept()