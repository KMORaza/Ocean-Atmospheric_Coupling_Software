from PyQt5.QtWidgets import QTextEdit
from PyQt5.QtGui import QFont

class ConsoleWidget(QTextEdit):
    def __init__(self):
        super().__init__()
        self.setReadOnly(True)
        self.setFont(QFont("Verdana", 10))
        self.setStyleSheet("background-color: #3c3f41; color: #ffffff;")
        self.append("Console initialized. Ready to run simulation.")
        
    def append_log(self, message):
        self.append(message)