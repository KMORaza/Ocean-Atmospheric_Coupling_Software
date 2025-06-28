import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget
from PyQt5.QtGui import QFont
from ControlPanel import ControlPanel
from PlotWidget import PlotWidget
from ConsoleWidget import ConsoleWidget

class MainApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Ocean-Atmosphere Coupling Model Simulator")
        self.setGeometry(100, 100, 1200, 800)
        
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
        """)
        
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)
        
        # Left side: Control panel and console
        left_layout = QVBoxLayout()
        self.control_panel = ControlPanel()
        self.console = ConsoleWidget()
        left_layout.addWidget(self.control_panel)
        left_layout.addWidget(self.console)
        
        # Right side: Plot widget
        self.plot_widget = PlotWidget()
        
        # Add layouts to main layout
        main_layout.addLayout(left_layout, 1)
        main_layout.addWidget(self.plot_widget, 2)
        
        # Connect control panel signals to slots
        self.control_panel.run_simulation.connect(self.run_simulation)
        
    def run_simulation(self, params):
        self.console.append_log("Starting simulation...")
        try:
            from Model import OceanAtmosphereModel
            model = OceanAtmosphereModel(**params)
            time, ocean_temps, atm_temps = model.run()
            self.plot_widget.update_plot(time, ocean_temps, atm_temps)
            self.console.append_log("Simulation completed successfully.")
        except Exception as e:
            self.console.append_log(f"Error: {str(e)}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setFont(QFont("Verdana", 10))
    window = MainApp()
    window.show()
    sys.exit(app.exec_())