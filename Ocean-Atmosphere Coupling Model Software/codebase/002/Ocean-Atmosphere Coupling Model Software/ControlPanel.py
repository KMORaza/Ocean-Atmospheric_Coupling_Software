from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QCheckBox
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtGui import QFont
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class ControlPanel(QWidget):
    start_simulation = pyqtSignal(dict)
    pause_simulation = pyqtSignal()
    reset_simulation = pyqtSignal()
    
    def __init__(self):
        super().__init__()
        logging.debug("Initializing ControlPanel")
        try:
            layout = QVBoxLayout()
            self.setLayout(layout)
            
            # Input fields with adjusted defaults for stability
            self.inputs = {}
            params = [
                ("Initial Ocean Temp (K)", "300.0"),
                ("Initial Atmosphere Temp (K)", "295.0"),
                ("Coupling Coefficient (W/m²/K)", "5.0"),  # Reduced for stability
                ("Ocean Heat Capacity (J/m²/K)", "4e6"),
                ("Atmosphere Heat Capacity (J/m²/K)", "1e5"),
                ("Time Step (s)", "1800.0"),  # Smaller time step
                ("Total Time (s)", "864000.0"),
                ("Grid Size", "50"),
                ("Coast Resolution Factor", "1.5"),  # Reduced for smoother transitions
                ("Nested Grid Size", "20"),
                ("AMR Threshold", "0.1"),  # Lower threshold for refinement
                ("Solar Forcing (W/m²)", "1361.0"),  # Realistic solar constant
                ("Longwave Coefficient", "0.6"),  # Adjusted for stability
                ("Advection Velocity (m/s)", "0.05")  # Reduced for stability
            ]
            
            for label, default in params:
                layout.addWidget(QLabel(label))
                edit = QLineEdit(default)
                edit.setFont(QFont("Verdana", 10))
                layout.addWidget(edit)
                self.inputs[label] = edit
            
            # Nested grid toggle
            self.use_nested_grid = QCheckBox("Enable Nested Grid")
            self.use_nested_grid.setFont(QFont("Verdana", 10))
            self.use_nested_grid.setChecked(False)
            layout.addWidget(self.use_nested_grid)
            
            # Control buttons
            start_button = QPushButton("Start Simulation")
            start_button.setFont(QFont("Verdana", 10))
            start_button.clicked.connect(self.emit_simulation_params)
            layout.addWidget(start_button)
            
            pause_button = QPushButton("Pause Simulation")
            pause_button.setFont(QFont("Verdana", 10))
            pause_button.clicked.connect(self.pause_simulation.emit)
            layout.addWidget(pause_button)
            
            reset_button = QPushButton("Reset Simulation")
            reset_button.setFont(QFont("Verdana", 10))
            reset_button.clicked.connect(self.reset_simulation.emit)
            layout.addWidget(reset_button)
            
            layout.addStretch()
            logging.debug("ControlPanel initialization complete")
        except Exception as e:
            logging.error(f"ControlPanel initialization failed: {str(e)}")
            raise
        
    def emit_simulation_params(self):
        logging.debug("Emitting simulation parameters")
        try:
            params = {
                "ocean_temp": float(self.inputs["Initial Ocean Temp (K)"].text()),
                "atm_temp": float(self.inputs["Initial Atmosphere Temp (K)"].text()),
                "coupling_coeff": float(self.inputs["Coupling Coefficient (W/m²/K)"].text()),
                "heat_capacity_ocean": float(self.inputs["Ocean Heat Capacity (J/m²/K)"].text()),
                "heat_capacity_atm": float(self.inputs["Atmosphere Heat Capacity (J/m²/K)"].text()),
                "time_step": float(self.inputs["Time Step (s)"].text()),
                "total_time": float(self.inputs["Total Time (s)"].text()),
                "grid_size": int(self.inputs["Grid Size"].text()),
                "coast_factor": float(self.inputs["Coast Resolution Factor"].text()),
                "use_nested_grid": self.use_nested_grid.isChecked(),
                "nested_grid_size": int(self.inputs["Nested Grid Size"].text()),
                "amr_threshold": float(self.inputs["AMR Threshold"].text()),
                "solar_forcing": float(self.inputs["Solar Forcing (W/m²)"].text()),
                "longwave_coeff": float(self.inputs["Longwave Coefficient"].text()),
                "adv_velocity": float(self.inputs["Advection Velocity (m/s)"].text())
            }
            # Validate inputs
            if params["grid_size"] < 10:
                raise ValueError("Grid Size must be at least 10.")
            if params["use_nested_grid"]:
                if params["nested_grid_size"] >= params["grid_size"]:
                    raise ValueError("Nested Grid Size must be smaller than Grid Size.")
                if params["nested_grid_size"] % 2 != 0:
                    raise ValueError("Nested Grid Size must be even.")
            if params["time_step"] <= 0 or params["total_time"] <= 0:
                raise ValueError("Time Step and Total Time must be positive.")
            if params["ocean_temp"] < 250 or params["ocean_temp"] > 350 or \
               params["atm_temp"] < 250 or params["atm_temp"] > 350:
                raise ValueError("Temperatures must be between 250K and 350K.")
            logging.debug("Simulation parameters validated")
            self.start_simulation.emit(params)
        except ValueError as e:
            logging.error(f"Input validation failed: {str(e)}")
            self.start_simulation.emit({"error": str(e)})