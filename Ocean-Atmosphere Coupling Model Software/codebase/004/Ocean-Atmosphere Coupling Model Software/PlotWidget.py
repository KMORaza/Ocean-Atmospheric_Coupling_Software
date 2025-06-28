import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QWidget, QVBoxLayout
from PyQt5.QtGui import QFont
import numpy as np
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class PlotWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        logging.debug("Initializing PlotWidget")
        try:
            self.figure = Figure(facecolor="#2b2b2b")
            self.canvas = FigureCanvas(self.figure)
            
            layout = QVBoxLayout()
            layout.addWidget(self.canvas)
            self.setLayout(layout)
            
            # Set font and dark theme
            font = {"family": "Verdana", "size": 10}
            matplotlib.rc("font", **font)
            matplotlib.rc("axes", edgecolor="white", labelcolor="white", facecolor="#3c3f41")
            matplotlib.rc("xtick", color="white")
            matplotlib.rc("ytick", color="white")
            matplotlib.rc("text", color="white")
            matplotlib.rc("figure", facecolor="#2b2b2b")
            
            # Initialize time series data for dynamic curve
            self.time_steps = []
            self.mean_ocean_temps = []
            self.mean_atm_temps = []
            
            logging.debug("PlotWidget initialization complete")
        except Exception as e:
            logging.error(f"PlotWidget initialization failed: {str(e)}")
            raise
    
    def update_plot(self, time, ocean_temps, atm_temps, refinement_mask, nested_offset, nested_size):
        logging.debug(f"Updating plot at time={time}")
        try:
            self.figure.clear()
            
            # Create subplots: 2 heatmaps (top) and 1 line plot (bottom)
            gs = self.figure.add_gridspec(2, 2, height_ratios=[3, 1])
            ax1 = self.figure.add_subplot(gs[0, 0])
            ax2 = self.figure.add_subplot(gs[0, 1])
            ax3 = self.figure.add_subplot(gs[1, :])
            
            # Plot heatmaps
            im1 = ax1.imshow(ocean_temps, cmap="viridis", origin="lower")
            im2 = ax2.imshow(atm_temps, cmap="plasma", origin="lower")
            
            # Overlay refinement mask
            ax1.contour(refinement_mask, colors="white", alpha=0.5, levels=[0.5])
            ax2.contour(refinement_mask, colors="white", alpha=0.5, levels=[0.5])
            
            # Overlay nested grid boundary
            if nested_offset is not None and nested_size is not None:
                rect = matplotlib.patches.Rectangle(
                    (nested_offset, nested_offset), nested_size, nested_size,
                    fill=False, edgecolor="yellow", linewidth=1
                )
                ax1.add_patch(rect)
                rect2 = matplotlib.patches.Rectangle(
                    (nested_offset, nested_offset), nested_size, nested_size,
                    fill=False, edgecolor="yellow", linewidth=1
                )
                ax2.add_patch(rect2)
            
            ax1.set_title(f"Ocean Temperature (K) at t={time:.2f}s")
            ax2.set_title(f"Atmosphere Temperature (K) at t={time:.2f}s")
            
            self.figure.colorbar(im1, ax=ax1)
            self.figure.colorbar(im2, ax=ax2)
            
            ax1.set_xlabel("X")
            ax1.set_ylabel("Y")
            ax2.set_xlabel("X")
            ax2.set_ylabel("Y")
            
            # Update time series data
            self.time_steps.append(time)
            self.mean_ocean_temps.append(np.mean(ocean_temps))
            self.mean_atm_temps.append(np.mean(atm_temps))
            
            # Plot dynamic curve
            ax3.plot(self.time_steps, self.mean_ocean_temps, label="Mean Ocean Temp (K)", color="cyan")
            ax3.plot(self.time_steps, self.mean_atm_temps, label="Mean Atm Temp (K)", color="magenta")
            ax3.set_title("Mean Temperatures Over Time")
            ax3.set_xlabel("Time (s)")
            ax3.set_ylabel("Temperature (K)")
            ax3.legend()
            ax3.grid(True, color="gray", linestyle="--", alpha=0.5)
            
            self.figure.tight_layout()
            self.canvas.draw()
            logging.debug("Plot update complete")
        except Exception as e:
            logging.error(f"Plot update failed: {str(e)}")
            raise
    
    def clear_plot(self):
        logging.debug("Clearing plot")
        try:
            self.figure.clear()
            self.time_steps = []
            self.mean_ocean_temps = []
            self.mean_atm_temps = []
            self.canvas.draw()
            logging.debug("Plot cleared")
        except Exception as e:
            logging.error(f"Plot clear failed: {str(e)}")
            raise