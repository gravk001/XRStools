from  PyQt4 import Qt, QtCore
import ui
import os
import PyMca5.PyMcaIO.EdfFile as edf

import PyMca5.PyMcaIO.specfile as specfile
import localfilesdialog
import myMaskImageWidget

import numpy

import roi
import spotdetection
import match

@ui.UILoadable
class spotdetectioncontrol(Qt.QDockWidget):
    def __init__(self, parent,flag,  detectionCallBack=None,  fatteningCallBack=None,
                 annotateMaskCallBack=None,    relabeliseMaskCallBack=None,
                 resetMask           = None):

        Qt.QDialog.__init__(self, parent, flag)
        self.loadUi()  # load the ui file

        self.detectionButton.clicked.connect(detectionCallBack)
        # self.fatteningSpinBox.valueChanged[int].connect(fatteningCallBack)

        self.fatteningCallBack = fatteningCallBack
        self.inflateButton.clicked.connect(self.callcallback)

        self.annotateButton.clicked.connect(annotateMaskCallBack)
        self.relabeliseButton.clicked.connect(relabeliseMaskCallBack)
        self.resetButton.clicked.connect(resetMask)



    def callcallback(self):
        value= self.fatteningSpinBox.value()
        self.fatteningCallBack(value)
   
@ui.UILoadable
class imageview(Qt.QWidget):
    def __init__(self, parent=None, isglobal=False):
        Qt.QWidget.__init__(self, parent)
        self.loadUi()  # load the ui file
        if not isglobal:
            self.registeringLayoutComboBox.addItems(["3x4 layout", "Vertical layout" ])

   
@ui.UILoadable
class spotregistrationcontrol(Qt.QDockWidget):
    def __init__(self, parent,flag,  globregistrationCallBack=None,
                 registrationCallBack=None):
        Qt.QDialog.__init__(self, parent, flag)
        self.loadUi()  # load the ui file
        
        self.globregistrationButton.clicked.connect(globregistrationCallBack)
        self.registrationButton.clicked.connect(registrationCallBack)

        self.tableWidget.setRowCount(3)
        self.tableWidget.setColumnCount(4)
        count=1
        for n in range(1,12+1):
            self.tableWidget.setItem( 2-(n-1)/4,(n-1)%4,  Qt.QTableWidgetItem("%d"%n))

        self.tableWidget.resizeColumnsToContents()

class MyTableModel(Qt.QAbstractTableModel):
    def __init__(self, parent, mylist, header, *args):
        Qt.QAbstractTableModel.__init__(self, parent, *args)
        self.mylist = mylist
        self.header = header
    def rowCount(self, parent):
        return len(self.mylist)
    def columnCount(self, parent):
        return len(self.mylist[0])
    def data(self, index, role):
        if not index.isValid():
            return None
        # elif role != Qt.DisplayRole:
        #     return None
        return self.mylist[index.row()][index.column()]
    def headerData(self, col, orientation, role):
        # if orientation == Qt.Horizontal and role == Qt.DisplayRole:
        #     return self.header[col]
        return None
    # def sort(self, col, order):
    #     """sort table by given column number col"""
    #     self.emit(SIGNAL("layoutAboutToBeChanged()"))
    #     self.mylist = sorted(self.mylist,
    #         key=operator.itemgetter(col))
    #     if order == Qt.DescendingOrder:
    #         self.mylist.reverse()
    #     self.emit(SIGNAL("layoutChanged()"))


@ui.UILoadable
class mainwindow(Qt.QMainWindow):
    def __init__(self, parent=None):
        Qt.QMainWindow.__init__(self, parent)
        self.loadUi()  # load the ui file


        self.actionLoad.triggered.connect(self.LoadRemote )
        self.actionSelectScanFiles.triggered.connect(self.LoadLocal)
        self.actionSpot_detection.triggered.connect(self.CreateSpotDetectionDockWidget)
        self.actionGlobalSpotDetection.triggered.connect(self.CreateGlobalSpotDetectionDockWidget)
        self.actionShowDatas.triggered.connect(self.showDatas)
        self.actionShowMasks.triggered.connect(self.showMasks)
        self.actionRegistration.triggered.connect(self.CreateRegistrationWidget)



 
    def detectionCallBack(self):
        print " in  detectionCallBack" 
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 

        globalMask = self.mws[0].getSelectionMask().astype("i")
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:itab-1], self.mws[1:itab]):
            offset += nofrois

        (name,geo,nofrois), mw   =     self.names_geos_nofrois[itab-1], self.mws[itab]
        self.detectSpotsSubMask( name,geo,nofrois, mw  , globalMask , offset )
        self.mws[0].setSelectionMask( globalMask )


    def annotateOneMaskCallBack(self):
        itab =  self.viewsTab.currentIndex()     
        if itab>0:
            self.mws[itab].annotateSpots(  )

    def annotateAllMasksCallBack(self):
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
            mw.annotateSpots(  )

    def relabeliseAllMasksCallBack(self):
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
            self.relabeliseSpots( mw, nofrois, name , geo, offset)

    def relabeliseOneMaskCallBack(self):
        itab =  self.viewsTab.currentIndex()     
        if itab==0: return
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:itab-1], self.mws[1:itab]):
            offset += nofrois
        (name,geo,nofrois), mw  =    self.names_geos_nofrois[itab-1], self.mws[itab]
        self.relabeliseSpots( mw , nofrois, name, geo, offset)

    def relabeliseSpots(self,mw, nofrois, name, geo, offset):

        mask = mw.getSelectionMask( ).astype("i")
        mask = (mask>0).astype("i")
        newmask = spotdetection.relabelise(mask,mask, nofrois)
        self.checkNspots(newmask.max(), nofrois , name) 
        mw.setSelectionMask( newmask )

        globalMask = self.mws[0].getSelectionMask().astype("i")
        globalMask[geo] = newmask
        self.mws[0].setSelectionMask( globalMask )

    def resetOneMask(self):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        mw = self.mws[itab]
        mw.graph.clearMarkers()
        mask = mw.getSelectionMask( ).astype("i")
        mask[:]=0
        mw.setSelectionMask( mask)

    def resetAllMasks(self):
        ret = self.warnForGloablChange()
        print ret
        if ret:
            for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
                mw.graph.clearMarkers()
                mask = mw.getSelectionMask( ).astype("i")
                mask[:]=0
                mw.setSelectionMask( mask)

    def fatteningCallBack(self,value):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        self.fatten(itab,value)

    def fatten(self, itab, value  ) :
        globalMask = self.mws[0].getSelectionMask().astype("i")
        name,geo,nofrois =   self.names_geos_nofrois[itab-1]

        mw = self.mws[itab]
        mask = mw.getSelectionMask( )
        mask = spotdetection.grow_mask(mask, value  ) 
        mw.setSelectionMask(mask  )

        globalMask[geo] = mask
        self.mws[0].setSelectionMask( globalMask )
        
    def checkNspots(self,nspots, nofrois,name ) :
        if nspots != nofrois:
            msgBox = Qt.QMessageBox ()
            msgBox.setText("Warning: found %d spots instead of expected %d "%(nspots , nofrois ) );
            msgBox.setInformativeText("For detector %s " %  name);
            msgBox.setStandardButtons(Qt.QMessageBox.Ok );
            msgBox.setDefaultButton(Qt.QMessageBox.Cancel);
            ret = msgBox.exec_();

    def    detectSpotsSubMask( self, name,geo,nofrois, mw   , globalMask, offset ):

        submatrix = mw.getImageData()
        
        mask  = spotdetection.get_spots_mask( submatrix, median_size = 5 ) 
    
        self.checkNspots(mask.max(), nofrois , name) 

        mw.setSelectionMask( mask)
        globalMask[geo]=mask+offset*numpy.less(0,mask)
        mw.annotateSpots(   )
        
    def GlobaldetectionCallBack(self, warn=False):

        if warn:
            print " in Global  detectionCallBack" 
            ret = self.warnForGloablChange()
        else:
            ret = True
        if ret:
            offset=0
            globalMask = self.mws[0].getSelectionMask().astype("i")
            for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois, self.mws[1:]):
                self.detectSpotsSubMask( name,geo,nofrois, mw  , globalMask , offset )
                offset += nofrois
            self.mws[0].setSelectionMask( globalMask )


    def showMasks(self):
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois, self.mws[1:]):
            mask = mw.getSelectionMask().astype("i")
            mw.setImageData(mask  , xScale=(0.0, 1.0), yScale=(0., 1.))

    def showDatas(self):
        Data = self.mws[0].getImageData()
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois, self.mws[1:]):
            d  = Data[geo]
            mw.setImageData(d  , xScale=(0.0, 1.0), yScale=(0., 1.))


    def warnForGloablChange(self):
        msgBox = Qt.QMessageBox ()
        msgBox.setText("You are going to recalculate the GLOBAL mask");
        msgBox.setInformativeText("This will reset all modifications to local masks. Do you want to proceed?");
        msgBox.setStandardButtons(Qt.QMessageBox.Ok |  Qt.QMessageBox.Cancel);
        msgBox.setDefaultButton(Qt.QMessageBox.Cancel);
        ret = msgBox.exec_();
        return ret==Qt.QMessageBox.Ok
        
    def GlobalfatteningCallBack(self,value):
        print " in  detectionCallBack ", value  
        ret = self.warnForGloablChange()
        print ret
        if ret:
            for itab in range(1, len(self.mws )  ) :
                self.fatten(itab,value)

    def CreateSpotDetectionDockWidget(self):
        w = spotdetectioncontrol(self, QtCore.Qt.Widget,detectionCallBack=self.detectionCallBack,  
                                 fatteningCallBack=self.fatteningCallBack,
                                 annotateMaskCallBack=self.annotateOneMaskCallBack,
                                 relabeliseMaskCallBack=self.relabeliseOneMaskCallBack,
                                 resetMask           = self.resetOneMask
                                 )
        self.addDockWidget ( QtCore.Qt.LeftDockWidgetArea, w )
        # w.setAllowedAreas (QtCore.Qt.AllDockWidgetAreas)
        w.show()

    def globregistration(self):
        print " in globregistration "
        for itab in range(1,len(self.mws)):
            self.registerTab(itab)

    def registration(self):
        itab =  self.viewsTab.currentIndex()     
        print self.layouts[itab].currentText()
        self.registerTab(itab)

    def registerTab(self, itab):
        if itab>0:
            name,geo,nofrois = self.names_geos_nofrois[itab-1]
            self.registerSpots( self.mws[itab], self.layouts[itab].currentText(), name, nofrois )
            self.mws[itab].annotateSpots( )

    def registerSpots( self,  mw, layoutString , name, nofrois) :

        mask = mw.getSelectionMask().astype("i")
        newmask = spotdetection.relabelise(mask,mask, nofrois)
        self.checkNspots(newmask.max(), nofrois , name) 
        mask = newmask
        nspots = mask.max()
        spots = []
        for i in range(1,nspots+1):
            zone = (mask==i)
            mask[zone]=mask[zone]+100
            m = zone.astype("f")
            msum=m.sum()
            if msum:
                ny,nx = m.shape
                px= (m.sum(axis=0)*numpy.arange(nx)).sum()/msum
                py= (m.sum(axis=1)*numpy.arange(ny)).sum()/msum
                spots.append((py,px,i))

        self.checkNspots(len(spots), nofrois , name) 
            
        print type(str(layoutString))
        if "3x4" in str(layoutString):
            
            positions = [ [px,py]   for (py,px,i) in spots ]
            choices = match.register(positions)
            newspots = []
            for (cx,cy),(y,x,i) in zip(choices, spots ):
                newspots.append((y,x,i,  (cy*4+cy)+1  ) )
        else:
            spots.sort()
            newspots = []
            for k,(y,x,i) in enumerate( spots):
                newspots.append((y,x,i,k+1) ) 

        for (y,x,i,k) in newspots:
            zone = (mask==(i+100))
            mask[zone] = k
                
        mw.setSelectionMask(mask)
            
    def CreateRegistrationWidget(self):

        w = spotregistrationcontrol(self, QtCore.Qt.Widget,globregistrationCallBack=self.globregistration,  
                                    registrationCallBack=self.registration
                                )

        self.addDockWidget ( QtCore.Qt.LeftDockWidgetArea, w )
        w.show()

    def CreateGlobalSpotDetectionDockWidget(self):
        w = spotdetectioncontrol(self, QtCore.Qt.Widget,detectionCallBack=self.GlobaldetectionCallBack, 
                                 fatteningCallBack=self.GlobalfatteningCallBack,
                                 annotateMaskCallBack=self.annotateAllMasksCallBack,
                                 relabeliseMaskCallBack=self.relabeliseAllMasksCallBack,
                                 resetMask           = self.resetAllMasks
                                 )
        w.detectionButton.setText("GLOBAL detection")
        self.addDockWidget ( QtCore.Qt.LeftDockWidgetArea, w )
        # w.setAllowedAreas (QtCore.Qt.AllDockWidgetAreas)
        w.show()

    def set_roiob(self, roiob):
        self.roiob=roiob

    def LoadRemote(self):
        if self.roiob is None:
            mb = qt.QMessageBox();
            mb.setText("No roi-object has been associated to the roi manager. Cannot Load.");
            mb.exec_();
            return
        self.showImage(self.roiob)


    def create_viewWidget(self, image = None, nofrois = None, changeTagOn=False, isglobal=False ) :
        view = imageview(self, isglobal)        
        maskW  =myMaskImageWidget.MaskImageWidget(self, aspect=True,
                                                  profileselection=True,
                                                  maxNRois=nofrois )
        maskW.setAcceptDrops(True)


        maskW.setDefaultColormap(2, logflag=True)

        maskW.changeTagOn = changeTagOn
        maskW.setImageData(image  , xScale=(0.0, 1.0), yScale=(0., 1.))
        maskW.setSelectionMask(image*0 , plot=True)
        view.roiContainerWidget.layout().addWidget(maskW)
        return view, maskW
            
    def showImage(self, roiob):

        subset_infos = roiob.get_geo_informations( roiob.image.shape )

        self.viewsTab.clear()       
        totNofrois = subset_infos["nofrois"]
        nofrois = totNofrois/len( subset_infos["subgeos"]  ) 

        self.mws=[]
        self.layouts=[]
        view, mw = self.create_viewWidget(image = roiob.image, nofrois = totNofrois , changeTagOn = False, isglobal=True  ) 
        self.viewsTab.addTab(view, "Global")
        self.mws.append(mw)
        self.layouts.append(view.registeringLayoutComboBox)

        self.names_geos_nofrois = zip(subset_infos["subnames"],  subset_infos["subgeos"], [nofrois]*len(subset_infos["subgeos"])  )

        for name,geo,nofr in self.names_geos_nofrois :
            view, mw = self.create_viewWidget(image = roiob.image[     geo ], nofrois = nofrois, changeTagOn = True  ) 
            self.viewsTab.addTab(view, name)
            self.mws.append(mw)
            self.layouts.append(view.registeringLayoutComboBox)
            

    def LoadLocal(self):
        if 0:
            w = localfilesdialog.localfilesdialog()
            w.exec_()
            sf =  str(w.SpecFileName_lineEdit.text())
            fn =  str(w.FileName_lineEdit.text())
            ns =  w.ScanNumber_spinBox.value()
        else:
            sf = "../data/test_files_scan425/raman"
            fn = "../data/test_files_scan425/h_raman_41125.edf"
            ns =  425

        template = getTemplateName( fn )

        s=specfile.Specfile(sf)
        Scan = s[ns-1]
        numbers = Scan.datacol("ccdno")

        roiob=None

        print numbers

        for n in numbers:

            image = [ edf.EdfFile(tok%n,"r").GetData(0)   for tok in template ] 
            image = numpy.concatenate(image, axis=0)
            if roiob is None:
                roiob = roi.rois( ) 
                shape = image.shape
                roiob.prepare_rois(  len(numbers) ,  "%dx%d"%shape ) 
            roiob.process(image)
    
        
        self.showImage(roiob)


        if 1:
            self.GlobaldetectionCallBack(warn=False)



def getTemplateName(name):
    dirname=os.path.dirname(str(name))
    name=os.path.basename(str(name))
    ls = len(name)
    fine = None
    inizio=None

    for n in range(ls)[::-1]:
        if fine is None and name[n] in "1234567890":
            fine = n+1
        if name[n] in  "1234567890":
            inizio=n
        if fine is not None and name[n] not in "1234567890":
            break
        
    print name
    print inizio
    print fine
        
    name = name[:inizio]+"%" + ("0%d"%(fine-inizio)+"d"+name[fine:])

    print name

    if name[:2] in ["h_", "v_"]:
        return [dirname+"/"+"h_"+name[2:],dirname+"/"+"v_"+name[2:]]
    else:
        return [dirname+"/"+name]


if __name__=="__main__":
    app=Qt.QApplication([])
    w = mainwindow()
    w.show()
    app.exec_()
