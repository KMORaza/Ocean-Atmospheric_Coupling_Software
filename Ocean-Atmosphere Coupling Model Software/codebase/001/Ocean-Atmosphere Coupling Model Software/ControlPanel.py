from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QCheckBox
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtGui import QFont

class ControlPanel(QWidget):
    run_simulation = pyqtSignal(dict)
    
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Input fields
        self.inputs = {}
        params = [
            ("Initial Ocean Temp (K)", "300.0"),
            ("Initial Atmosphere Temp (K)", "295.0"),
            ("Coupling Coefficient (W/m²/K)", "10.0"),
            ("Ocean Heat Capacity (J/m²/K)", "4e6"),
            ("Atmosphere Heat Capacity (J/m²/K)", "1e5"),
            ("regular: ("Time Step (s)", "3600.0"),
            ("Total Time (s)", "864000.0"),
            ("Grid Size", "50"),
            ("Coast Resolution Factor", "2.0"),
            ("Nested Grid Size", "20"),
            ("AMR Threshold", "0.5")
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
        
        # Run button
        run_button = QPushButton("Run Simulation")
        run_button.setFont(QFont("Verdana", 10))
        run_button.clicked.connect(self.emit_simulation_params)
        layout.addWidget(run_button)
        
        layout.addStretch()
        
    def emit_simulation_params(self):
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
                "amr_threshold": float(self.inputs["AMR Threshold"].text())
            }
            self.run_simulation.emit(params)
        except ValueError as e:
            # Error handling will be displayed in console
            pass