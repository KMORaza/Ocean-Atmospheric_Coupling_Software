import numpy as np
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGroupBox, QFormLayout, QWidget, QMessageBox, QSizePolicy
from PyQt5.QtGui import QFont, QPainter, QImage, QPen, QBrush, QColor
from PyQt5.QtCore import Qt, QTimer
import logging

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

    def initialize_field(self, grid_size, eddy_strength, eddy_radius, num_eddies):
        self.grid_size = grid_size
        x = np.linspace(-1, 1, grid_size)
        y = np.linspace(-1, 1, grid_size)
        X, Y = np.meshgrid(x, y)
        self.vorticity = np.zeros((grid_size, grid_size))
        
        # Initialize multiple eddies with random positions
        self.eddy_centers = []
        self.eddy_strengths = []
        self.eddy_radii = []
        np.random.seed(42)  # For reproducibility
        for _ in range(num_eddies):
            x0 = np.random.uniform(-0.8, 0.8)
            y0 = np.random.uniform(-0.8, 0.8)
            strength = eddy_strength * np.random.uniform(0.8, 1.2)
            radius = eddy_radius * np.random.uniform(0.8, 1.2)
            self.eddy_centers.append([x0, y0])
            self.eddy_strengths.append(strength)
            self.eddy_radii.append(radius)
            # Rankine vortex model
            r = np.sqrt((X - x0)**2 + (Y - y0)**2)
            vortex = np.where(r < radius, strength * (r / radius), strength * (radius / r))
            self.vorticity += vortex * np.exp(-r**2 / (2 * radius**2))
        
        self.update_visualization()

    def update_simulation(self, vorticity_diffusion, rotation_rate, background_flow, decay_rate):
        if self.vorticity is None:
            return
        dt = 0.1
        self.time += dt
        x = np.linspace(-1, 1, self.grid_size)
        X, Y = np.meshgrid(x, x)
        
        # Update eddy positions (translation due to background flow)
        for center in self.eddy_centers:
            center[0] += background_flow * dt  # Move eddies with background flow
            if center[0] > 1:  # Periodic boundary
                center[0] -= 2
        
        # Recalculate vorticity field with updated positions and decay
        self.vorticity *= np.exp(-decay_rate * dt)  # Exponential decay
        new_vorticity = np.zeros_like(self.vorticity)
        for (x0, y0), strength, radius in zip(self.eddy_centers, self.eddy_strengths, self.eddy_radii):
            r = np.sqrt((X - x0)**2 + (Y - y0)**2)
            vortex = np.where(r < radius, strength * (r / radius), strength * (radius / r))
            new_vorticity += vortex * np.exp(-r**2 / (2 * radius**2))
        
        # Add diffusion
        self.vorticity = new_vorticity + vorticity_diffusion * dt * (
            np.roll(new_vorticity, 1, axis=0) + 
            np.roll(new_vorticity, -1, axis=0) +
            np.roll(new_vorticity, 1, axis=1) +
            np.roll(new_vorticity, -1, axis=1) - 
            4 * new_vorticity
        )
        
        # Add background shear flow
        self.vorticity += background_flow * (Y * 0.1)  # Linear shear
        self.update_visualization()

    def update_visualization(self):
        self.image.fill(Qt.white)
        painter = QPainter(self.image)
        if self.vorticity is not None:
            norm_vorticity = (self.vorticity - np.min(self.vorticity)) / (np.max(self.vorticity) - np.min(self.vorticity) + 1e-10)
            for i in range(self.grid_size):
                for j in range(self.grid_size):
                    value = norm_vorticity[i, j]
                    # Smooth color gradient from cyan to blue
                    r = int(0 * value)
                    g = int(255 * (1 - value))
                    b = int(255 * value)
                    color = QColor(r, g, b)
                    painter.setPen(QPen(color))
                    painter.setBrush(QBrush(color))
                    x = int(400 * i / self.grid_size)
                    y = int(400 * j / self.grid_size)
                    painter.drawRect(x, y, 4, 4)
            
            # Draw streamlines (simplified)
            painter.setPen(QPen(Qt.black, 1))
            for _ in range(10):  # Draw 10 streamlines
                x0 = np.random.uniform(-1, 1)
                y0 = np.random.uniform(-1, 1)
                points = []
                for _ in range(20):
                    i = int((x0 + 1) * self.grid_size / 2)
                    j = int((y0 + 1) * self.grid_size / 2)
                    if 0 <= i < self.grid_size and 0 <= j < self.grid_size:
                        vx = self.vorticity[j, i] * 0.1  # Approximate velocity
                        vy = -self.vorticity[j, i] * 0.1
                        points.append((int(400 * (x0 + 1) / 2), int(400 * (y0 + 1) / 2)))
                        x0 += vx * 0.01
                        y0 += vy * 0.01
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
            self.setGeometry(200, 200, 600, 800)
            self.setFont(QFont("Consolas", 10))
            
            # Apply modern mobile-like stylesheet with 90s aesthetic
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
            """)
            
            self.setup_ui()
            self.timer = QTimer()
            self.timer.timeout.connect(self.update_simulation)
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(False)
            logging.debug("OceanicEddyAndFrontWindow initialized")
        except Exception as e:
            logging.error(f"OceanicEddyAndFrontWindow initialization failed: {str(e)}")
            raise

    def setup_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)
        self.setLayout(main_layout)
        
        # Simulation widget (takes most space)
        self.simulation_widget = EddySimulationWidget()
        self.simulation_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(self.simulation_widget, stretch=3)
        
        # Control panel (card-like, at bottom)
        control_group = QGroupBox("Eddy Simulation Parameters")
        control_group.setStyleSheet("QGroupBox { background-color: #FFFFFF; box-shadow: 0px 2px 5px rgba(0,0,0,0.2); }")
        form_layout = QFormLayout()
        form_layout.setHorizontalSpacing(10)
        form_layout.setVerticalSpacing(8)
        
        self.inputs = {}
        params = [
            ("Grid Size", "100", "Grid resolution (10–200)"),
            ("Eddy Strength", "1.0", "Vortex intensity (0.5–2.0)"),
            ("Eddy Radius", "0.2", "Vortex size (0.1–0.3)"),
            ("Number of Eddies", "3", "Number of vortices (1–10)"),
            ("Vorticity Diffusion", "0.01", "Diffusion rate (0.005–0.02)"),
            ("Rotation Rate (rad/s)", "0.1", "Spin speed (0.05–0.2)"),
            ("Background Flow (m/s)", "0.05", "Flow speed (0.02–0.1)"),
            ("Decay Rate (1/s)", "0.01", "Dissipation rate (0.005–0.02)")
        ]
        
        for label, default, tooltip in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Consolas", 10))
            edit.setMinimumWidth(100)
            edit.setToolTip(tooltip)
            form_layout.addRow(QLabel(label), edit)
            self.inputs[label] = edit
        
        control_group.setLayout(form_layout)
        main_layout.addWidget(control_group, stretch=1)
        
        # Buttons (in a horizontal layout, mobile-like)
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
            eddy_strength = float(self.inputs["Eddy Strength"].text())
            eddy_radius = float(self.inputs["Eddy Radius"].text())
            num_eddies = int(self.inputs["Number of Eddies"].text())
            
            if grid_size < 10:
                raise ValueError("Grid Size must be at least 10")
            if eddy_strength <= 0:
                raise ValueError("Eddy Strength must be positive")
            if eddy_radius <= 0 or eddy_radius > 1:
                raise ValueError("Eddy Radius must be between 0 and 1")
            if num_eddies < 1 or num_eddies > 10:
                raise ValueError("Number of Eddies must be between 1 and 10")
                
            self.simulation_widget.initialize_field(grid_size, eddy_strength, eddy_radius, num_eddies)
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            logging.debug("Eddy simulation initialized")
            
        except ValueError as e:
            logging.error(f"Initialization failed: {str(e)}")
            QMessageBox.warning(self, "Input Error", str(e))

    def start_simulation(self):
        try:
            if self.simulation_widget.vorticity is None:
                self.initialize_simulation()
            self.timer.start(100)
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(True)
            self.reset_button.setEnabled(True)
            logging.debug("Eddy simulation started")
        except Exception as e:
            logging.error(f"Start simulation failed: {str(e)}")
            QMessageBox.warning(self, "Start Error", str(e))

    def pause_simulation(self):
        try:
            self.timer.stop()
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            logging.debug("Eddy simulation paused")
        except Exception as e:
            logging.error(f"Pause simulation failed: {str(e)}")
            QMessageBox.warning(self, "Pause Error", str(e))

    def reset_simulation(self):
        try:
            self.timer.stop()
            self.simulation_widget.vorticity = None
            self.simulation_widget.time = 0.0
            self.simulation_widget.eddy_centers = []
            self.simulation_widget.eddy_strengths = []
            self.simulation_widget.eddy_radii = []
            self.simulation_widget.image.fill(Qt.white)
            self.simulation_widget.update()
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.initialize_simulation()
            logging.debug("Eddy simulation reset")
        except Exception as e:
            logging.error(f"Reset failed: {str(e)}")
            QMessageBox.warning(self, "Reset Error", str(e))

    def update_simulation(self):
        try:
            vorticity_diffusion = float(self.inputs["Vorticity Diffusion"].text())
            rotation_rate = float(self.inputs["Rotation Rate (rad/s)"].text())
            background_flow = float(self.inputs["Background Flow (m/s)"].text())
            decay_rate = float(self.inputs["Decay Rate (1/s)"].text())
            self.simulation_widget.update_simulation(vorticity_diffusion, rotation_rate, background_flow, decay_rate)
        except ValueError as e:
            logging.error(f"Simulation update failed: {str(e)}")
            QMessageBox.warning(self, "Simulation Error", str(e))
            self.pause_simulation()

    def closeEvent(self, event):
        self.timer.stop()
        logging.debug("OceanicEddyAndFrontWindow closed")
        event.accept()
