import sys
import logging
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QMessageBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QTimer
try:
    from ControlPanel import ControlPanel
    from PlotWidget import PlotWidget
    from ConsoleWidget import ConsoleWidget
    from Model import OceanAtmosphereModel
except ImportError as e:
    print(f"Import error: {str(e)}")
    raise

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class MainApp(QMainWindow):
    def __init__(self):
        super().__init__()
        logging.debug("Initializing MainApp")
        try:
            self.setWindowTitle("Ocean-Atmosphere Coupling Model Simulator")
            self.setGeometry(100, 100, 1200, 900)
            
            # Set default font
            self.setFont(QFont("Verdana", 10))
            
            # Apply dark theme stylesheet
            self.setStyleSheet("""
                QMainWindow, QWidget { background-color: #2b2b2b; color: #ffffff; }
                QLineEdit { background-color: #3c3f41; color: #ffffff; border: 1px solid #555555; }
                QPushButton { background-color: #4a4a4a; color: #ffffff; border: 1px solid #555555; }
                QPushButton:hover { background-color: #5a5a5a; }
                QTextEdit { background-color: #3c3f41; color: #ffffff; }
                QLabel { color: #ffffff; }
                QCheckBox { color: #ffffff; }
            """)
            
            # Main widget and layout
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_layout = QHBoxLayout()
            main_widget.setLayout(main_layout)
            
            # Left side: Control panel and console
            left_layout = QVBoxLayout()
            logging.debug("Creating ControlPanel")
            self.control_panel = ControlPanel()
            logging.debug("Creating ConsoleWidget")
            self.console = ConsoleWidget()
            left_layout.addWidget(self.control_panel)
            left_layout.addWidget(self.console)
            
            # Right side: Plot widget
            logging.debug("Creating PlotWidget")
            self.plot_widget = PlotWidget()
            
            # Add layouts to main layout
            main_layout.addLayout(left_layout, 1)
            main_layout.addWidget(self.plot_widget, 2)
            
            # Simulation state
            self.model = None
            self.timer = QTimer()
            self.timer.timeout.connect(self.update_simulation)
            self.current_step = 0
            
            # Connect control panel signals
            self.control_panel.start_simulation.connect(self.start_simulation)
            self.control_panel.pause_simulation.connect(self.pause_simulation)
            self.control_panel.reset_simulation.connect(self.reset_simulation)
            
            logging.debug("MainApp initialization complete")
        except Exception as e:
            logging.error(f"MainApp initialization failed: {str(e)}")
            QMessageBox.critical(None, "Initialization Error", f"Failed to start application: {str(e)}")
            raise
    
    def start_simulation(self, params):
        if "error" in params:
            self.console.append_log(f"Input error: {params['error']}")
            return
        self.console.append_log("Starting simulation...")
        try:
            # Validate nested grid size
            if params["use_nested_grid"]:
                if params["nested_grid_size"] >= params["grid_size"]:
                    raise ValueError("Nested grid size must be smaller than global grid size.")
                if params["nested_grid_size"] % 2 != 0:
                    raise ValueError("Nested grid size must be even for 2x refinement.")
            self.model = OceanAtmosphereModel(**params)
            self.current_step = 0
            self.timer.start(100)  # Update every 100ms
            self.console.append_log("Simulation started.")
        except Exception as e:
            self.console.append_log(f"Error: {str(e)}")
            logging.error(f"Start simulation failed: {str(e)}")
            self.model = None
        
    def pause_simulation(self):
        self.timer.stop()
        self.console.append_log("Simulation paused.")
        
    def reset_simulation(self):
        self.timer.stop()
        self.model = None
        self.current_step = 0
        self.plot_widget.clear_plot()
        self.console.append_log("Simulation reset.")
        
    def update_simulation(self):
        if self.model and self.current_step < int(self.model.total_time / self.model.dt):
            try:
                time, ocean_temps, atm_temps, refinement_mask = self.model.step(self.current_step)
                self.plot_widget.update_plot(time, ocean_temps, atm_temps, refinement_mask,
                                          self.model.nested_grid.offset if self.model.nested_grid else None,
                                          self.model.nested_grid.nested_grid_size if self.model.nested_grid else None)
                self.console.append_log(f"Step {self.current_step}: Time = {time:.2f} s, "
                                      f"Mean Ocean Temp = {np.mean(ocean_temps):.2f} K, "
                                      f"Mean Atm Temp = {np.mean(atm_temps):.2f} K")
                self.current_step += 1
            except Exception as e:
                self.timer.stop()
                self.console.append_log(f"Simulation error at step {self.current_step}: {str(e)}")
                logging.error(f"Simulation step failed: {str(e)}")
        else:
            self.timer.stop()
            self.console.append_log("Simulation completed.")

if __name__ == "__main__":
    logging.debug("Starting application")
    try:
        app = QApplication(sys.argv)
        app.setFont(QFont("Verdana", 10))
        window = MainApp()
        window.show()
        logging.debug("Application window shown")
        sys.exit(app.exec_())
    except Exception as e:
        logging.error(f"Application failed to start: {str(e)}")
        print(f"Application failed to start: {str(e)}")
        raise