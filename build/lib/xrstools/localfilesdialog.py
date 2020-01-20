from  PyQt4 import Qt
import ui



import PyMca5.PyMcaIO.specfile as specfile




# this will use the ui module to load
# the structure specified in localfilesdialog.ui
@ui.UILoadable
class localfilesdialog(Qt.QDialog):
    def __init__(self, parent=None):
        Qt.QDialog.__init__(self, parent)
        self.loadUi()  # load the ui file
        self.BrowseSpec_pushButton.clicked.connect(self.__onBrowseSpec)
        self.BrowseImage_pushButton.clicked.connect(self.__onBrowseFile)
        self.SpecFileName_lineEdit.textChanged.connect(self.__onChangeSpec)
        self.ScanNumber_spinBox.setMaximum(-1)

    def __onBrowseSpec(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        self.SpecFileName_lineEdit.setText(filename)

    def __onChangeSpec(self):
        filename =  str(self.SpecFileName_lineEdit.text())
        print filename

        try:
            s=specfile.Specfile(filename)
        except:
            s=None
        if s is not None:
            ns = len(s)
            self.ScanNumber_spinBox.setMinimum(0)
            self.ScanNumber_spinBox.setMaximum(ns-1)
        else:
            self.ScanNumber_spinBox.setMaximum(-1)

# s[450].alllabels()
# s[450]["ccdno"]
# s[424].datacol("ccdno")
# s[425].datacol("ccdno")
# ls
# history




   
    def __onBrowseFile(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        self.FileName_lineEdit.setText(filename)
   



if __name__=="__main__":
    app=Qt.QApplication([])
    w = localfilesdialog()
    w.show()
    app.exec_()
