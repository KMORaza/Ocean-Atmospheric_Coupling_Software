from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, 
                             QPushButton, QCheckBox, QScrollArea, QGroupBox, 
                             QFormLayout, QMessageBox, QSizePolicy)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QFont, QPalette, QColor
import logging
import traceback

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class ControlPanel(QWidget):
    start_simulation = pyqtSignal(dict)
    pause_simulation = pyqtSignal()
    reset_simulation = pyqtSignal()
    set_time_scales = pyqtSignal(float, float)
    open_turbulent_mixing = pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        logging.debug("Initializing ControlPanel")
        try:
            # Main layout
            main_layout = QVBoxLayout()
            main_layout.setContentsMargins(10, 10, 10, 10)
            main_layout.setSpacing(15)
            self.setLayout(main_layout)
            
            # Create scroll area
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QScrollArea.NoFrame)
            scroll_content = QWidget()
            scroll_layout = QVBoxLayout(scroll_content)
            scroll_layout.setContentsMargins(5, 5, 5, 5)
            scroll_layout.setSpacing(15)
            
            # Group parameters into logical sections
            self.create_initial_conditions_group(scroll_layout)
            self.create_physical_parameters_group(scroll_layout)
            self.create_time_scales_group(scroll_layout)
            self.create_numerical_parameters_group(scroll_layout)
            self.create_grid_settings_group(scroll_layout)
            self.create_flux_parameters_group(scroll_layout)
            self.create_conservation_parameters_group(scroll_layout)
            
            scroll.setWidget(scroll_content)
            scroll.setMaximumHeight(500)
            main_layout.addWidget(scroll)
            
            self.create_control_buttons(main_layout)
            
            # Connect pause button to update button states
            self.pause_simulation.connect(self.on_pause_simulation)
            
            logging.debug("ControlPanel initialization complete")
        except Exception as e:
            logging.error(f"ControlPanel initialization failed: {str(e)}\n{traceback.format_exc()}")
            QMessageBox.critical(self, "Initialization Error", 
                               f"Failed to initialize control panel: {str(e)}")
            raise
        
    def create_initial_conditions_group(self, layout):
        group = QGroupBox("Initial Conditions")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        self.inputs = {}
        params = [
            ("Initial Ocean Temp (K)", "300.0"),
            ("Initial Atmosphere Temp (K)", "295.0"),
            ("Coupling Coefficient (W/m²/K)", "5.0")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            edit.setMinimumWidth(150)
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_physical_parameters_group(self, layout):
        group = QGroupBox("Physical Parameters")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        params = [
            ("Ocean Heat Capacity (J/m²/K)", "4e6"),
            ("Atmosphere Heat Capacity (J/m²/K)", "1e5"),
            ("Solar Forcing (W/m²)", "1361.0"),
            ("Longwave Coefficient", "0.6"),
            ("Advection Velocity (m/s)", "0.05"),
            ("Drag Coefficient", "0.0012"),
            ("Wind Speed (m/s)", "5.0")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_time_scales_group(self, layout):
        group = QGroupBox("Time Scale Coupling")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        params = [
            ("Ocean Time Scale Factor (slower)", "1.0"),
            ("Atmosphere Time Scale Factor (faster)", "0.1")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            tooltip = ("Larger values make ocean evolve slower\n"
                      "Typical range: 0.5-2.0 for ocean, 0.05-0.2 for atmosphere")
            edit.setToolTip(tooltip)
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_numerical_parameters_group(self, layout):
        group = QGroupBox("Numerical Parameters")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        params = [
            ("Time Step (s)", "1800.0"),
            ("Total Time (s)", "864000.0")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_grid_settings_group(self, layout):
        group = QGroupBox("Grid Settings")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        params = [
            ("Grid Size", "50"),
            ("Coast Resolution Factor", "1.5"),
            ("Nested Grid Size", "20"),
            ("AMR Threshold", "0.1")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        self.use_nested_grid = QCheckBox("Enable Nested Grid")
        self.use_nested_grid.setFont(QFont("Verdana", 10))
        form.addRow(self.use_nested_grid)
        
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_flux_parameters_group(self, layout):
        group = QGroupBox("Flux Parameters")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        params = [
            ("Precipitation Rate (kg/m²/s)", "1e-6"),
            ("Evaporation Rate (kg/m²/s)", "1e-6"),
            ("Turbulent Mixing Coefficient (m²/s)", "1e-4"),
            ("CO2 Transfer Velocity (m/s)", "1e-5")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_conservation_parameters_group(self, layout):
        group = QGroupBox("Conservation Parameters")
        form = QFormLayout()
        form.setHorizontalSpacing(15)
        
        params = [
            ("Freshwater Conservation Coeff", "1.0"),
            ("CO2 Conservation Coeff", "1.0")
        ]
        
        for label, default in params:
            edit = QLineEdit(default)
            edit.setFont(QFont("Verdana", 10))
            form.addRow(QLabel(label), edit)
            self.inputs[label] = edit
            
        group.setLayout(form)
        layout.addWidget(group)
    
    def create_control_buttons(self, layout):
        button_layout = QHBoxLayout()
        button_layout.setContentsMargins(0, 15, 0, 0)
        button_layout.setSpacing(10)
        
        self.start_button = QPushButton("Start Simulation")
        self.start_button.setFont(QFont("Verdana", 10, QFont.Bold))
        self.start_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.start_button.setMinimumHeight(40)
        self.start_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border-radius: 5px;
                padding: 8px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:pressed {
                background-color: #3d8b40;
            }
            QPushButton:disabled {
                background-color: #a5d6a7;
            }
        """)
        self.start_button.clicked.connect(self.emit_simulation_params)
        
        self.pause_button = QPushButton("Pause")
        self.pause_button.setFont(QFont("Verdana", 10, QFont.Bold))
        self.pause_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.pause_button.setMinimumHeight(40)
        self.pause_button.setStyleSheet("""
            QPushButton {
                background-color: #FFC107;
                color: black;
                border-radius: 5px;
                padding: 8px;
            }
            QPushButton:hover {
                background-color: #e6ac00;
            }
            QPushButton:pressed {
                background-color: #cc9900;
            }
            QPushButton:disabled {
                background-color: #ffe082;
            }
        """)
        self.pause_button.clicked.connect(self.pause_simulation.emit)
        self.pause_button.setEnabled(False)
        
        self.reset_button = QPushButton("↻ Reset")
        self.reset_button.setFont(QFont("Verdana", 10, QFont.Bold))
        self.reset_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.reset_button.setMinimumHeight(40)
        self.reset_button.setStyleSheet("""
            QPushButton {
                background-color: #F44336;
                color: white;
                border-radius: 5px;
                padding: 8px;
            }
            QPushButton:hover {
                background-color: #e53935;
            }
            QPushButton:pressed {
                background-color: #d32f2f;
            }
            QPushButton:disabled {
                background-color: #ef9a9a;
            }
        """)
        self.reset_button.clicked.connect(self.reset_simulation.emit)
        
        self.turbulent_mixing_button = QPushButton("Turbulent Mixing")
        self.turbulent_mixing_button.setFont(QFont("Verdana", 10, QFont.Bold))
        self.turbulent_mixing_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.turbulent_mixing_button.setMinimumHeight(40)
        self.turbulent_mixing_button.setStyleSheet("""
            QPushButton {
                background-color: #0288D1;
                color: white;
                border-radius: 5px;
                padding: 8px;
            }
            QPushButton:hover {
                background-color: #0277BD;
            }
            QPushButton:pressed {
                background-color: #01579B;
            }
            QPushButton:disabled {
                background-color: #4FC3F7;
            }
        """)
        self.turbulent_mixing_button.clicked.connect(self.open_turbulent_mixing.emit)
        self.turbulent_mixing_button.setEnabled(False)
        
        button_layout.addWidget(self.start_button)
        button_layout.addWidget(self.pause_button)
        button_layout.addWidget(self.reset_button)
        button_layout.addWidget(self.turbulent_mixing_button)
        layout.addLayout(button_layout)
    
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
                "adv_velocity": float(self.inputs["Advection Velocity (m/s)"].text()),
                "drag_coeff": float(self.inputs["Drag Coefficient"].text()),
                "wind_speed": float(self.inputs["Wind Speed (m/s)"].text()),
                "precip_rate": float(self.inputs["Precipitation Rate (kg/m²/s)"].text()),
                "evap_rate": float(self.inputs["Evaporation Rate (kg/m²/s)"].text()),
                "mixing_coeff": float(self.inputs["Turbulent Mixing Coefficient (m²/s)"].text()),
                "co2_transfer_coeff": float(self.inputs["CO2 Transfer Velocity (m/s)"].text()),
                "freshwater_conservation_coeff": float(self.inputs["Freshwater Conservation Coeff"].text()),
                "co2_conservation_coeff": float(self.inputs["CO2 Conservation Coeff"].text()),
                "ocean_time_scale": float(self.inputs["Ocean Time Scale Factor (slower)"].text()),
                "atm_time_scale": float(self.inputs["Atmosphere Time Scale Factor (faster)"].text())
            }
            
            self.validate_parameters(params)
            self.set_time_scales.emit(params["ocean_time_scale"], params["atm_time_scale"])
            
            # Update button states
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(True)
            self.reset_button.setEnabled(True)
            self.turbulent_mixing_button.setEnabled(True)
            
            logging.debug("Simulation parameters validated")
            self.start_simulation.emit(params)
            
        except ValueError as e:
            error_msg = f"Invalid input: {str(e)}"
            logging.error(error_msg)
            QMessageBox.warning(self, "Input Error", error_msg)
            self.start_simulation.emit({"error": str(e)})
        except Exception as e:
            error_msg = f"Unexpected error: {str(e)}\n{traceback.format_exc()}"
            logging.error(error_msg)
            QMessageBox.critical(self, "Error", error_msg)
            self.start_simulation.emit({"error": str(e)})
    
    def on_pause_simulation(self):
        logging.debug("Pause simulation triggered")
        try:
            # Update button states to allow restarting
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.turbulent_mixing_button.setEnabled(True)
        except Exception as e:
            logging.error(f"Pause simulation handling failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to handle pause: {str(e)}")
    
    def validate_parameters(self, params):
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
        if params["drag_coeff"] <= 0 or params["wind_speed"] < 0 or \
           params["precip_rate"] < 0 or params["evap_rate"] < 0 or \
           params["mixing_coeff"] <= 0 or params["co2_transfer_coeff"] <= 0 or \
           params["freshwater_conservation_coeff"] <= 0 or params["co2_conservation_coeff"] <= 0:
            raise ValueError("Physical coefficients must be non-negative.")
        if params["ocean_time_scale"] <= 0 or params["atm_time_scale"] <= 0:
            raise ValueError("Time scale factors must be positive.")
        if params["ocean_time_scale"] <= params["atm_time_scale"]:
            raise ValueError("Ocean time scale should be larger (slower) than atmosphere time scale.")