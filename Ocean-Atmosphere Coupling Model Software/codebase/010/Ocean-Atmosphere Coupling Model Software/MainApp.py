import sys
import logging
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QMessageBox, QScrollArea
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QTimer

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Import modules with error handling
try:
    from ControlPanel import ControlPanel
except ImportError as e:
    logging.error(f"Failed to import ControlPanel: {str(e)}")
    raise
try:
    from PlotWidget import PlotWidget
except ImportError as e:
    logging.error(f"Failed to import PlotWidget: {str(e)}")
    raise
try:
    from ConsoleWidget import ConsoleWidget
except ImportError as e:
    logging.error(f"Failed to import ConsoleWidget: {str(e)}")
    raise
try:
    from Model import OceanAtmosphereModel
except ImportError as e:
    logging.error(f"Failed to import Model: {str(e)}")
    raise
try:
    from TurbulentMixing import TurbulentMixingWindow
except ImportError as e:
    logging.error(f"Failed to import TurbulentMixing: {str(e)}")
    raise
try:
    from SurfaceLayerPhysics import SurfaceLayerPhysicsWindow
except ImportError as e:
    logging.error(f"Failed to import SurfaceLayerPhysics: {str(e)}")
    raise
try:
    from BoundaryLayerSchemes import BoundaryLayerSchemesWindow
except ImportError as e:
    logging.error(f"Failed to import BoundaryLayerSchemes: {str(e)}")
    raise
try:
    from HighOrderTimeStepping import HighOrderTimeSteppingWindow
except ImportError as e:
    logging.error(f"Failed to import HighOrderTimeStepping: {str(e)}")
    raise
try:
    from AirSeaInteraction import AirSeaInteractionWindow
except ImportError as e:
    logging.error(f"Failed to import AirSeaInteraction: {str(e)}")
    raise
try:
    from CloudMicroPhysics import CloudMicroPhysicsWindow
except ImportError as e:
    logging.error(f"Failed to import CloudMicroPhysics: {str(e)}")
    print(f"Warning: CloudMicroPhysics module not found. Cloud Microphysics functionality will be disabled: {str(e)}")
    CloudMicroPhysicsWindow = None
try:
    from OceanicEddyAndFront import OceanicEddyAndFrontWindow
except ImportError as e:
    logging.error(f"Failed to import OceanicEddyAndFront: {str(e)}")
    print(f"Warning: OceanicEddyAndFront module not found. Oceanic Eddies functionality will be disabled: {str(e)}")
    OceanicEddyAndFrontWindow = None

class MainApp(QMainWindow):
    def __init__(self):
        super().__init__()
        logging.debug("Initializing MainApp")
        try:
            self.setWindowTitle("Ocean-Atmosphere Coupling Model Simulator")
            self.setGeometry(100, 100, 1200, 900)
            
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
                QTextEdit { 
                    background-color: #FFFFFF; 
                    color: #000000; 
                    border: 2px inset #808080; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QLabel { 
                    color: #000000; 
                    padding: 2px; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QCheckBox { 
                    color: #000000; 
                    font-family: Consolas;
                    font-size: 10pt;
                }
                QScrollArea { 
                    background-color: #C0C0C0; 
                    border: 2px inset #808080; 
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
            
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_layout = QHBoxLayout()
            main_widget.setLayout(main_layout)
            
            left_layout = QVBoxLayout()
            logging.debug("Creating ControlPanel")
            self.control_panel = ControlPanel()
            logging.debug("Creating ConsoleWidget")
            self.console = ConsoleWidget()
            
            console_scroll = QScrollArea()
            console_scroll.setWidgetResizable(True)
            console_scroll.setWidget(self.console)
            console_scroll.setMinimumHeight(200)
            console_scroll.setMaximumHeight(1000)
            
            left_layout.addWidget(self.control_panel)
            left_layout.addWidget(console_scroll)
            left_layout.addStretch()
            
            logging.debug("Creating PlotWidget")
            self.plot_widget = PlotWidget()
            
            main_layout.addLayout(left_layout, 1)
            main_layout.addWidget(self.plot_widget, 2)
            
            self.model = None
            self.timer = QTimer()
            self.timer.timeout.connect(self.update_simulation)
            self.current_step = 0
            self.turbulent_mixing_window = None
            self.surface_layer_physics_window = None
            self.boundary_layer_schemes_window = None
            self.high_order_time_stepping_window = None
            self.air_sea_interaction_window = None
            self.cloud_microphysics_window = None
            self.oceanic_eddies_window = None
            
            self.control_panel.start_simulation.connect(self.start_simulation)
            self.control_panel.pause_simulation.connect(self.pause_simulation)
            self.control_panel.reset_simulation.connect(self.reset_simulation)
            self.control_panel.open_turbulent_mixing.connect(self.open_turbulent_mixing)
            self.control_panel.open_surface_layer_physics.connect(self.open_surface_layer_physics)
            self.control_panel.open_boundary_layer_schemes.connect(self.open_boundary_layer_schemes)
            self.control_panel.open_high_order_time_stepping.connect(self.open_high_order_time_stepping)
            self.control_panel.air_sea_interaction.connect(self.open_air_sea_interaction)
            if CloudMicroPhysicsWindow is not None:
                self.control_panel.cloud_microphysics.connect(self.open_cloud_microphysics)
            else:
                self.control_panel.cloud_microphysics_button.setEnabled(False)
                self.console.append_log("Warning: Cloud Microphysics module not found. Functionality disabled.")
            if OceanicEddyAndFrontWindow is not None:
                self.control_panel.oceanic_eddies.connect(self.open_oceanic_eddies)
                self.control_panel.oceanic_eddies_button.setEnabled(True)  # MODIFIED: Explicitly enable button
                logging.debug("Oceanic Eddies signal connected")
            else:
                self.control_panel.oceanic_eddies_button.setEnabled(False)
                self.console.append_log("Warning: Oceanic Eddies module not found. Functionality disabled.")
            
            logging.debug("MainApp initialization complete")
        except Exception as e:
            logging.error(f"MainApp initialization failed: {str(e)}")
            QMessageBox.critical(self, "Initialization Error", f"Failed to start application: {str(e)}")
            raise
    
    def start_simulation(self, params):
        if "error" in params:
            self.console.append_log(f"Input error: {params['error']}")
            return
        self.console.append_log("Starting simulation...")
        try:
            if params["use_nested_grid"]:
                if params["nested_grid_size"] >= params["grid_size"]:
                    raise ValueError("Nested grid size must be smaller than global grid size.")
                if params["nested_grid_size"] % 2 != 0:
                    raise ValueError("Nested grid size must be even for 2x refinement.")
            self.model = OceanAtmosphereModel(**params)
            self.current_step = 0
            self.timer.start(100)
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
        if self.turbulent_mixing_window is not None:
            self.turbulent_mixing_window.close()
            self.turbulent_mixing_window = None
            self.control_panel.turbulent_mixing_button.setEnabled(False)
        if self.surface_layer_physics_window is not None:
            self.surface_layer_physics_window.close()
            self.surface_layer_physics_window = None
            self.control_panel.slp_button.setEnabled(False)
        if self.boundary_layer_schemes_window is not None:
            self.boundary_layer_schemes_window.close()
            self.boundary_layer_schemes_window = None
            self.control_panel.bls_button.setEnabled(False)
        if self.high_order_time_stepping_window is not None:
            self.high_order_time_stepping_window.close()
            self.high_order_time_stepping_window = None
            self.control_panel.hots_button.setEnabled(False)
        if self.air_sea_interaction_window is not None:
            self.air_sea_interaction_window.close()
            self.air_sea_interaction_window = None
            self.control_panel.air_sea_button.setEnabled(True)
        if self.cloud_microphysics_window is not None:
            self.cloud_microphysics_window.close()
            self.cloud_microphysics_window = None
            self.control_panel.cloud_microphysics_button.setEnabled(CloudMicroPhysicsWindow is not None)
        if self.oceanic_eddies_window is not None:
            self.oceanic_eddies_window.close()
            self.oceanic_eddies_window = None
            self.control_panel.oceanic_eddies_button.setEnabled(OceanicEddyAndFrontWindow is not None)
        
    def open_turbulent_mixing(self):
        logging.debug("Opening turbulent mixing window")
        try:
            if self.model is None:
                self.console.append_log("Error: No active simulation. Start simulation first.")
                QMessageBox.warning(self, "Error", "Please start a simulation before opening the turbulent mixing window.")
                return
            if self.turbulent_mixing_window is None:
                self.turbulent_mixing_window = TurbulentMixingWindow(self.model, parent=self)
            self.turbulent_mixing_window.show()
            self.console.append_log("Turbulent mixing window opened.")
        except Exception as e:
            self.console.append_log(f"Error opening turbulent mixing window: {str(e)}")
            logging.error(f"Open turbulent mixing failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to open turbulent mixing window: {str(e)}")
    
    def open_surface_layer_physics(self):
        logging.debug("Opening surface layer physics window")
        try:
            if self.model is None:
                self.console.append_log("Error: No active simulation. Start simulation first.")
                QMessageBox.warning(self, "Error", "Please start a simulation before opening the surface layer physics window.")
                return
            if self.surface_layer_physics_window is None:
                self.surface_layer_physics_window = SurfaceLayerPhysicsWindow(self.model, parent=self)
            self.surface_layer_physics_window.show()
            self.console.append_log("Surface layer physics window opened.")
        except Exception as e:
            self.console.append_log(f"Error opening surface layer physics window: {str(e)}")
            logging.error(f"Open surface layer physics failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to open surface layer physics window: {str(e)}")
    
    def open_boundary_layer_schemes(self):
        logging.debug("Opening boundary layer schemes window")
        try:
            if self.model is None:
                self.console.append_log("Error: No active simulation. Start simulation first.")
                QMessageBox.warning(self, "Error", "Please start a simulation before opening the boundary layer schemes window.")
                return
            if self.boundary_layer_schemes_window is None:
                self.boundary_layer_schemes_window = BoundaryLayerSchemesWindow(self.model, parent=self)
            self.boundary_layer_schemes_window.show()
            self.console.append_log("Boundary layer schemes window opened.")
        except Exception as e:
            self.console.append_log(f"Error opening boundary layer schemes window: {str(e)}")
            logging.error(f"Open boundary layer schemes failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to open boundary layer schemes window: {str(e)}")
    
    def open_high_order_time_stepping(self):
        logging.debug("Opening high-order time stepping window")
        try:
            if self.model is None:
                self.console.append_log("Error: No active simulation. Start simulation first.")
                QMessageBox.warning(self, "Error", "Please start a simulation before opening the high-order time stepping window.")
                return
            if self.high_order_time_stepping_window is None:
                self.high_order_time_stepping_window = HighOrderTimeSteppingWindow(self.model, parent=self)
            self.high_order_time_stepping_window.show()
            self.console.append_log("High-order time stepping window opened.")
        except Exception as e:
            self.console.append_log(f"Error opening high-order time stepping window: {str(e)}")
            logging.error(f"Open high-order time stepping failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to open high-order time stepping window: {str(e)}")
    
    def open_air_sea_interaction(self):
        logging.debug("Opening air-sea interaction window")
        try:
            if self.air_sea_interaction_window is None:
                self.air_sea_interaction_window = AirSeaInteractionWindow(parent=self)
            self.air_sea_interaction_window.show()
            self.console.append_log("Air-sea interaction window opened.")
        except Exception as e:
            self.console.append_log(f"Error opening air-sea interaction window: {str(e)}")
            logging.error(f"Open air-sea interaction failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to open air-sea interaction window: {str(e)}")
    
    def open_cloud_microphysics(self):
        logging.debug("Opening cloud microphysics window")
        try:
            if CloudMicroPhysicsWindow is None:
                self.console.append_log("Error: Cloud Microphysics module not available.")
                QMessageBox.warning(self, "Error", "Cloud Microphysics module not found. Please ensure CloudMicroPhysics.py is in the project directory.")
                return
            if self.cloud_microphysics_window is None:
                self.cloud_microphysics_window = CloudMicroPhysicsWindow(parent=self)
            self.cloud_microphysics_window.show()
            self.console.append_log("Cloud microphysics window opened.")
        except Exception as e:
            self.console.append_log(f"Error opening cloud microphysics window: {str(e)}")
            logging.error(f"Open cloud microphysics failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to open cloud microphysics window: {str(e)}")
    
    def open_oceanic_eddies(self):
        logging.debug("Attempting to open oceanic eddies window")  # MODIFIED: Already present, kept for clarity
        print("Attempting to open oceanic eddies window")  # MODIFIED: Added for immediate console feedback
        try:
            if OceanicEddyAndFrontWindow is None:
                self.console.append_log("Error: Oceanic Eddies module not available.")
                QMessageBox.warning(self, "Error", "Oceanic Eddies module not found. Please ensure OceanicEddyAndFront.py is in the project directory.")
                logging.warning("Oceanic Eddies module not found")
                return
            logging.debug("Creating OceanicEddyAndFrontWindow instance")  # MODIFIED: Added debug log
            print("Creating OceanicEddyAndFrontWindow instance")  # MODIFIED: Added console output
            if self.oceanic_eddies_window is None or not self.oceanic_eddies_window.isVisible():
                self.oceanic_eddies_window = OceanicEddyAndFrontWindow(parent=self)
                logging.debug("Showing oceanic eddies window")  # MODIFIED: Added debug log
                print("Showing oceanic eddies window")  # MODIFIED: Added console output
                self.oceanic_eddies_window.show()
                self.oceanic_eddies_window.move(self.frameGeometry().center())  # MODIFIED: Center the window
                self.console.append_log("Oceanic eddies window opened.")
                logging.debug("Oceanic eddies window opened successfully")
            else:
                self.oceanic_eddies_window.raise_()
                self.oceanic_eddies_window.activateWindow()
                self.console.append_log("Oceanic eddies window brought to front.")
                logging.debug("Oceanic eddies window raised")
        except Exception as e:
            error_msg = f"Error opening oceanic eddies window: {str(e)}"
            self.console.append_log(error_msg)
            logging.error(error_msg)
            print(error_msg)  # MODIFIED: Added console output
            QMessageBox.critical(self, "Error", f"Failed to open oceanic eddies window: {str(e)}")
    
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

    def closeEvent(self, event):
        logging.debug("Closing MainApp")
        try:
            self.timer.stop()
            self.model = None
            if self.turbulent_mixing_window:
                self.turbulent_mixing_window.close()
            if self.surface_layer_physics_window:
                self.surface_layer_physics_window.close()
            if self.boundary_layer_schemes_window:
                self.boundary_layer_schemes_window.close()
            if self.high_order_time_stepping_window:
                self.high_order_time_stepping_window.close()
            if self.air_sea_interaction_window:
                self.air_sea_interaction_window.close()
            if self.cloud_microphysics_window:
                self.cloud_microphysics_window.close()
            if self.oceanic_eddies_window:
                self.oceanic_eddies_window.close()
            self.console.append_log("Application closed.")
            logging.debug("MainApp closed successfully")
            event.accept()
        except Exception as e:
            logging.error(f"Close MainApp failed: {str(e)}")
            self.console.append_log(f"Error: Failed to close application: {str(e)}")
            event.accept()

if __name__ == "__main__":
    logging.debug("Starting application")
    try:
        app = QApplication(sys.argv)
        app.setFont(QFont("Consolas", 10))
        window = MainApp()
        window.show()
        logging.debug("Application window shown")
        sys.exit(app.exec_())
    except Exception as e:
        logging.error(f"Application failed to start: {str(e)}")
        print(f"Application failed to start: {str(e)}")
        raise