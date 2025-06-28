from PyQt5.QtWidgets import QTextEdit
from PyQt5.QtGui import QFont
import logging

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class ConsoleWidget(QTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        logging.debug("Initializing ConsoleWidget")
        try:
            self.setReadOnly(True)
            self.setFont(QFont("Verdana", 10))
            self.setStyleSheet("background-color: #3c3f41; color: #ffffff;")
            logging.debug("ConsoleWidget initialization complete")
        except Exception as e:
            logging.error(f"ConsoleWidget initialization failed: {str(e)}")
            raise
    
    def append_log(self, message):
        logging.debug(f"Appending log: {message}")
        try:
            self.append(message)
        except Exception as e:
            logging.error(f"Console append failed: {str(e)}")
            raise