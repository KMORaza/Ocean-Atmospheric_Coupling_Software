import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QWidget, QVBoxLayout
from PyQt5.QtGui import QFont
import numpy as np

class PlotWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(facecolor="#2b2b2b")
        self.canvas = FigureCanvas(self.figure)
        
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        
        # Set font and dark theme for plot
        font = {"family": "Verdana", "size": 10, "color": "white"}
        matplotlib.rc("font", **font)
        matplotlib.rc("axes", edgecolor="white", labelcolor="white", facecolor="#3c3f41")
        matplotlib.rc("xtick", color="white")
        matplotlib.rc("ytick", color="white")
        matplotlib.rc("text", color="white")
        
    def update_plot(self, time, ocean_temps, atm_temps):
        self.figure.clear()
        
        # Create subplots for ocean and atmosphere heatmaps
        ax1 = self.figure.add_subplot(121)
        ax2 = self.figure.add_subplot(122)
        
        # Plot final time step as heatmap
        im1 = ax1.imshow(ocean_temps[-1], cmap="viridis", origin="lower")
        im2 = ax2.imshow(atm_temps[-1], cmap="plasma", origin="lower")
        
        ax1.set_title("Ocean Temperature (K)")
        ax2.set_title("Atmosphere Temperature (K)")
        
        self.figure.colorbar(im1, ax=ax1)
        self.figure.colorbar(im2, ax=ax2)
        
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        
        self.figure.tight_layout()
        self.canvas.draw()