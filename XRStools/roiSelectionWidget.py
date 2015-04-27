from  PyQt4 import Qt, QtCore
import ui
import os
import PyMca5.PyMcaIO.EdfFile as edf
import PyMca5

import PyMca5.PyMcaIO.specfile as specfile
import localfilesdialog
import myMaskImageWidget

import numpy
import string

try:
    import xrs_rois
except:
    import roi as xrs_rois

import spotdetection
import match
import pickle
import h5py


FASTDEBUG=False
try:
    import PyTango
except:
    print " Could not load PyTango"


@ui.UILoadable
class spotdetectioncontrol(Qt.QDockWidget):
    def __init__(self, parent,flag,  detectionCallBack=None,  fatteningCallBack=None,
                 thresholdCallBack=None,
                 annotateMaskCallBack=None,    relabeliseMaskCallBack=None,
                 resetMask           = None):

        Qt.QDialog.__init__(self, parent, flag)
        self.loadUi()  # load the ui file

        self.detectionButton.clicked.connect(detectionCallBack)
        # self.fatteningSpinBox.valueChanged[int].connect(fatteningCallBack)

        self.fatteningCallBack = fatteningCallBack
        self.inflateButton.clicked.connect(self.callcallback)


        self.thresholdCallBack = thresholdCallBack
        self.thresholdButton.clicked.connect(self.callThrecallback)


        self.annotateButton.clicked.connect(annotateMaskCallBack)
        self.relabeliseButton.clicked.connect(relabeliseMaskCallBack)
        self.resetButton.clicked.connect(resetMask)

        self.geo_informations = None

    def callThrecallback(self):
        print self.thresholdEdit.text()
        print str(self.thresholdEdit.text())
        value= string.atof(str(self.thresholdEdit.text()))
        self.thresholdCallBack(value)

    def callcallback(self):
        value= self.fatteningSpinBox.value()
        self.fatteningCallBack(value)
   
@ui.UILoadable
class imageview(Qt.QWidget):
    all_layouts=[]
    def __init__(self, parent=None, isglobal=False):
        Qt.QWidget.__init__(self, parent)
        self.loadUi()  # load the ui file
        if not isglobal:
            self.registeringLayoutComboBox.addItems(["3x4 layout", "Vertical layout" ])
            self.all_layouts.append(self.registeringLayoutComboBox)
        else:
            self.registeringLayoutComboBox.addItems(["3x4 layout", "Vertical layout" ])
            self.registeringLayoutComboBox.currentIndexChanged.connect(self.changeAll)
    
    def changeAll(self,index):
        for t in self.all_layouts:
            t.setCurrentIndex(index)

   
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

        def _mousePressEvent( event):
            """ in case of problems, generate the code for spotregistrationcontrol from ui
                and reimplement it as a class method """

            print " evento at ", event.pos()
            print "  widget ",   self.tableWidget.itemAt( event.pos())
            print    self.tableWidget.itemAt( event.pos()).text()


            drag      =  Qt.QDrag(   self.tableWidget)
            mimeData  =  Qt.QMimeData()

            mimeData.setText( self.tableWidget.itemAt( event.pos()).text()   );
            drag.setMimeData(mimeData)
            ## iconPixmap = Qt.QPixmap(_cross_data)
            ## drag.setPixmap(iconPixmap)

            dropAction = drag.exec_();


        self.tableWidget.mousePressEvent = _mousePressEvent
        

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
    def __init__(self, parent=None, labelformat="ROI%02d"):


        self.geo_informations=None
        self.isOK=0
        self.labelformat = labelformat

        Qt.QMainWindow.__init__(self, parent )
        self.loadUi()  # load the ui file


        # self.actionLoad.triggered.connect(self.LoadRemote )
        self.actionSelectScanFiles.triggered.connect(self.LoadLocal)
        self.actionSpot_detection.triggered.connect(self.CreateSpotDetectionDockWidget)
        self.actionGlobalSpotDetection.triggered.connect(self.CreateGlobalSpotDetectionDockWidget)
        self.actionShowDatas.triggered.connect(self.showDatas)
        self.actionShowMasks.triggered.connect(self.showMasks)
        self.actionRegistration.triggered.connect(self.CreateRegistrationWidget)


        self.actionWrite_mask_on_file.triggered.connect(  self.write_mask_on_file  )
        self.actionLoad_mask_from_file.triggered.connect(  self.read_mask_from_file  )
        ## self.actionRemote_load.triggered.connect(  self.remoteMaskload  )
        self.actionRemote_load.triggered.connect(  self.remoteMaskload  )
        self.actionPush_mask_remotely.triggered.connect(  self.PushMask  )
        self.actionWrite_masksDict_on_file.triggered.connect(  self.write_masksDict_on_file  )
        self.actionLoad_masksDict_from_file.triggered.connect(  self.load_masksDict_from_file  )
        
        
        self.actionExit.triggered.connect(self.Exit)
        self.actionConfirm_And_Exit.triggered.connect(self.confirm_and_Exit)

        self.actionLoad_image_from_hdf5.triggered.connect(self.load_image_from_hdf5)

        self.actionLoad_maskDict_from_hdf5.triggered.connect(self.load_maskDict_from_hdf5)




    def Exit(self):
        self.isOK=0
        self.close()
    def confirm_and_Exit(self):
        self.isOK=1
        self.close()

    def remoteMaskload(self):
        ## 
        pathtolima = "id20/xrs/mpx-ram"
        dev = PyTango.DeviceProxy(pathtolima ) 

        string =dev.getRoifrompickle()
        masks = pickle.loads(string)
        masksDict = self.masks2MasksDict(masks)

        self.load_masksDict(masksDict)
            

    def PushMask(self):
        ## 
        pathtolima = "id20/xrs/mpx-ram"
        dev = PyTango.DeviceProxy(pathtolima ) 
 
        masks = self.getMasks()

        print " de mainWindow je vais pusher : " , masks
        print " ========================================" 
        masks_string = pickle.dumps(masks)
        open("mask","w").write(masks_string)
        dev.setRoifromPickle( masks_string  )

    def getMasks(self) :
        totnofrois = self.recomposeGlobalMask()
        globalMask = self.mws[0].getSelectionMask().astype("i")
        masks = [] 
        nrois = globalMask.max()
        nummaxroi = (globalMask).max()
        masks = []
        for n in range(1,totnofrois +1):
            rawindices = numpy.nonzero( globalMask == n )
            if(len(rawindices)):
                if(len(rawindices[0])):
                    Y1 = rawindices[0].min()
                    Y2 = rawindices[0].max()+1
                    X1 = rawindices[1].min()
                    X2 = rawindices[1].max()+1
                    submask = (globalMask[Y1:Y2, X1:X2 ] == n).astype("i")
                    # submask = (globalMask[Y1:Y2, X1:X2 ] / n).astype("i"
                    masks.append(  ( n, [Y1,X1], submask  )  )
            else:
                masks.append(  ( n, [-1,-1], numpy.array([[1]])  )  )
        return masks

    # def getMasks(self) :
    #     totnofrois = self.recomposeGlobalMask()
    #     globalMask = self.mws[0].getSelectionMask().astype("i")
    #     masks = [] 
    #     nrois = globalMask.max()
    #     nummaxroi = (globalMask).max()
    #     masks = []
    #     for n in range(1, nummaxroi+1):
    #         rawindices = numpy.nonzero( globalMask == n )
    #         if(len(rawindices)):
    #             if(len(rawindices[0])):
    #                 Y1 = rawindices[0].min()
    #                 Y2 = rawindices[0].max()+1
    #                 X1 = rawindices[1].min()
    #                 X2 = rawindices[1].max()+1
    #                 submask = (globalMask[Y1:Y2, X1:X2 ] == n).astype("i")
    #                 masks.append(  ( n, [Y1,X1], submask  )  )
    #     return masks

    def getMasksDict(self) :
        masks = self.getMasks()
        return self.masks2MasksDict(masks)

    def masks2MasksDict(self,masks):
        masksDict={}
        for m in masks:
            if(m[1][0]>=0 and m[1][1]>=0  ):
                masksDict[ self.labelformat% (m[0]-1)  ]=[m[1],m[2]]
        return masksDict
        


    def recomposeGlobalMask(self):
      offset=0
      globalMask = self.mws[0].getSelectionMask().astype("i")
      for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
          localmask = numpy.less(0,mw.getSelectionMask() ) 
          globalMask[geo] =  offset*localmask +  mw.getSelectionMask()
          offset += nofrois
      self.mws[0].setSelectionMask(globalMask  )
      return offset

    
    def decomposeGlobalMask(self):
      offset=0
      globalMask = self.mws[0].getSelectionMask().astype("i")
      for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
          localmask = numpy.less(0,globalMask[geo]  ) 
          Mask =    globalMask[geo] - offset *localmask
          offset += nofrois
          mw.setSelectionMask(Mask  )
          

    def write_mask_on_file(self):
        filename =  Qt.QFileDialog.getSaveFileName()
        print filename
        if filename is not None:
            filename=str(filename)
            self.recomposeGlobalMask()
            globalMask = self.mws[0].getSelectionMask().astype("i")
            ef = edf.EdfFile( filename, "w+")
            ef.WriteImage( {}, globalMask )

    def write_masksDict_on_file(self):
        filename =  Qt.QFileDialog.getSaveFileName()
        print filename
        if filename is not None:
            filename=str(filename)
            masksDict = self.getMasksDict()
            filename=str(filename)
            f = open(filename, 'wb')
            pickle.dump(masksDict ,  f)
            f.close()
            
    def load_masksDict_from_file(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        print filename
        if filename is not None:
            filename=str(filename)
            f = open(filename, 'rb')
            masksDict = pickle.load( f)
            f.close()
            self.load_masksDict( masksDict)


    def load_masksDict(self, masksDict):
            self.recomposeGlobalMask()
            mask  = self.mws[0].getSelectionMask().astype("i")
            mask[:]=0
            mask =   convert_redmatrix_to_matrix(masksDict,mask, offsetX=0, offsetY=0)
            self.mws[0].setSelectionMask(mask)
            self.decomposeGlobalMask()

    def read_mask_from_file(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        if filename is not None:
            filename=str(filename)
            ef = edf.EdfFile( filename, "r")
            mask = ef.GetData(0)
            self.recomposeGlobalMask()
            self.mws[0].setSelectionMask(mask)
            self.decomposeGlobalMask()
 
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

    def get_geo(self):
        if self.geo_informations is None:
            subset_infos = xrs_rois.get_geo_informations( self.image.shape )
        else:
            subset_infos = self.geo_informations 
        return subset_infos

    def annotateOneMaskCallBack(self):
        itab =  self.viewsTab.currentIndex()     
        if itab>0:
            subset_infos = self.get_geo()

            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

            self.mws[itab].annotateSpots( a_ids )

    def annotateAllMasksCallBack(self):
        itab=1
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):

            subset_infos = self.get_geo()
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

            mw.annotateSpots( a_ids )
            itab+=1

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



    def threshold(self,itab,value ):
        globalMask = self.mws[0].getSelectionMask().astype("i")
        name,geo,nofrois =   self.names_geos_nofrois[itab-1]

        mw = self.mws[itab]
        mask = mw.getSelectionMask( )
        print mask.sum()

        data = self.image[geo]

        mask = spotdetection.threshold_mask(mask, data ,  value  ) 
        mw.setSelectionMask(mask  )

        globalMask[geo] = mask
        self.mws[0].setSelectionMask( globalMask )
        
        print mask.sum()


    def localThresholdCallBack(self,value):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        self.threshold(itab,value)

    def globalThresholdCallBack(self,value):
        ret = self.warnForGloablChange()
        print ret
        if ret:
            for itab in range(1, len(self.mws )  ) :
                self.threshold(itab,value)


    def fatten(self, itab, value  ) :
        globalMask = self.mws[0].getSelectionMask().astype("i")
        name,geo,nofrois =   self.names_geos_nofrois[itab-1]

        mw = self.mws[itab]
        mask = mw.getSelectionMask( )
        mask = spotdetection.grow_mask(mask, value  ) 
        mw.setSelectionMask(mask  )

        globalMask[geo] = mask
        self.mws[0].setSelectionMask( globalMask )


    def fatteningCallBack(self,value):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        self.fatten(itab,value)

    def GlobalfatteningCallBack(self,value):
        print " in  detectionCallBack ", value  
        ret = self.warnForGloablChange()
        print ret
        if ret:
            for itab in range(1, len(self.mws )  ) :
                self.fatten(itab,value)


        
    def checkNspots(self,nspots, nofrois,name ) :
        if nspots != nofrois:
            msgBox = Qt.QMessageBox ()
            msgBox.setText("Warning: found %d spots instead of expected %d "%(nspots , nofrois ) );
            msgBox.setInformativeText("For detector %s " %  name);
            msgBox.setStandardButtons(Qt.QMessageBox.Ok );
            msgBox.setDefaultButton(Qt.QMessageBox.Cancel);
            ret = msgBox.exec_();

    def    detectSpotsSubMask( self, name,geo,nofrois, mw   , globalMask, offset ):



        itab =  self.viewsTab.currentIndex()
        if itab: 
            subset_infos = self.get_geo() 
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

                
        else:
            a_ids = None

        submatrix = mw.getImageData()
        
        mask  = spotdetection.get_spots_mask( submatrix, median_size = 5 ) 
    
        self.checkNspots(mask.max(), nofrois , name) 

        mw.setSelectionMask( mask)
        globalMask[geo]=mask+offset*numpy.less(0,mask)



        mw.annotateSpots( a_ids   )
        
    def GlobaldetectionCallBack(self, warn=True):
        print " GlobaldetectionCallBack ",  warn
        if 1 or warn:
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
        

    def CreateSpotDetectionDockWidget(self):
        w = spotdetectioncontrol(self, QtCore.Qt.Widget,detectionCallBack=self.detectionCallBack,  
                                 fatteningCallBack=self.fatteningCallBack,
                                 thresholdCallBack  =   self.localThresholdCallBack, 
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

            subset_infos = self.get_geo() 
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]]["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

            self.mws[itab].annotateSpots(a_ids )

    def registerSpots( self,  mw, layoutString , name, nofrois) :

        mask = mw.getSelectionMask().astype("i")
        newmask = spotdetection.relabelise(mask,mask, nofrois)
        self.checkNspots(newmask.max(), nofrois , name) 
        mask = newmask
        nspots = mask.max()
        spots = []
        for i in range(1,nspots+1):
            zone = (mask==i)
            mask[zone]=i+100
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
            print str( numpy.array(positions))
            choices = match.register(numpy.array(positions)  )
            

            newspots = []
            for (cx,cy),(y,x,i) in zip(choices, spots ):
                print x,y, cx, cy
                newspots.append((y,x,i,  (cy*4+cx)+1  ) )
        else:
            spots.sort()
            newspots = []
            for k,(y,x,i) in enumerate( spots):
                newspots.append((y,x,i,k+1) ) 

        print " NEWSPOTS ", newspots

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
                                 thresholdCallBack  =   self.globalThresholdCallBack, 
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
            
    def showImage(self, image, geo_informations=None):
        self.image=image
        if geo_informations is None:
            geo_informations = self.get_geo()  
        self.geo_informations = geo_informations

        self.image=image

        subset_infos = geo_informations

        self.viewsTab.clear()       
        totNofrois = subset_infos["nofrois"]
        nofrois = totNofrois/len( subset_infos["subgeos"]  ) 

        self.mws=[]
        self.layouts=[]
        view, mw = self.create_viewWidget(image = image, nofrois = totNofrois , changeTagOn = False, isglobal=True  ) 
        self.viewsTab.addTab(view, "Global")
        self.mws.append(mw)
        self.layouts.append(view.registeringLayoutComboBox)

        self.names_geos_nofrois = zip(subset_infos["subnames"],  subset_infos["subgeos"], [nofrois]*len(subset_infos["subgeos"])  )

        for name,geo,nofr in self.names_geos_nofrois :
            view, mw = self.create_viewWidget(image = image[     geo ], nofrois = nofrois, changeTagOn = True  ) 
            self.viewsTab.addTab(view, name)
            self.mws.append(mw)
            self.layouts.append(view.registeringLayoutComboBox)
            

    def LoadLocal(self):
        if not FASTDEBUG:
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
                roiob = xrs_rois.rois( ) 
                shape = image.shape
                roiob.prepare_rois(  len(numbers) ,  "%dx%d"%shape ) 
            roiob.process(image)
    
        self.geo_informations=None
        self.showImage(roiob.image)
        self.roiob = roiob

        if FASTDEBUG:
            self.GlobaldetectionCallBack(warn=False)


    def load_maskDict_from_hdf5(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        filename=str(filename)
        print filename
        if len(filename):
            import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget
            storage=[None]
            def mySlot(ddict):
                name = ddict["name"]
                storage[0]=name
            # browse
            self.__hdf5Dialog = hdf5dialog()
            self.__hdf5Dialog.setWindowTitle('Select your data set by a double click')
            self.__hdf5Dialog.mainLayout =  self.__hdf5Dialog.verticalLayout_2
            fileModel = HDF5Widget.FileModel()
            fileView = HDF5Widget.HDF5Widget(fileModel)
            hdf5File = fileModel.openFile(filename)
            shiftsDataset = None
            fileView.sigHDF5WidgetSignal.connect(mySlot)
            self.__hdf5Dialog.mainLayout.addWidget(fileView)
            self.__hdf5Dialog.resize(400, 200)

            ret = self.__hdf5Dialog.exec_()
            print ret
            hdf5File.close()
            
            if ret:
                print " Obtained " 
                name =  storage[0]
                print name
                file= h5py.File(filename,"r")
                datagroup = file[name]
                masks={}
                xrs_rois.load_rois_fromh5(datagroup,masks)
                file.close()
                self.load_masksDict(masks)

    def load_image_from_hdf5(self):
        print " load " 
        filename =  Qt.QFileDialog.getOpenFileName()
        print " OK " 
        filename=str(filename)
        print filename
        if len(filename):
            import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget

            storage=[None]

            def mySlot(ddict):
                name = ddict["name"]
                storage[0]=name
                print " MY SLOT " 
                print name

            # browse
            self.__hdf5Dialog = hdf5dialog()
            self.__hdf5Dialog.setWindowTitle('Select your data set by a double click')
            self.__hdf5Dialog.mainLayout =  self.__hdf5Dialog.verticalLayout_2
            fileModel = HDF5Widget.FileModel()
            fileView = HDF5Widget.HDF5Widget(fileModel)
            hdf5File = fileModel.openFile(filename)
            shiftsDataset = None
            fileView.sigHDF5WidgetSignal.connect(mySlot)
            self.__hdf5Dialog.mainLayout.addWidget(fileView)
            self.__hdf5Dialog.resize(400, 200)
            
            # self.__hdf5Dialog.setModal(True)
            # self.__hdf5Dialog.show()

            ret = self.__hdf5Dialog.exec_()

            hdf5File.close()
            
            if ret:
                print " Obtained " 
                name =  storage[0]
                file= h5py.File(filename,"r")
                image4roi = file[name][:]
                file.close()

                self.showImage(image4roi)



    

@ui.UILoadable
class hdf5dialog(Qt.QDialog):
    def __init__(self, parent=None):
        Qt.QDialog.__init__(self, parent)
        self.loadUi()  # load the ui file
       

def convert_redmatrix_to_matrix( masksDict,mask, offsetX=0, offsetY=0):
    for key, (pos,M)  in masksDict.iteritems():
        num=int("".join([c for c in key if c.isdigit()]))
        S = M.shape
        inset =    (slice(offsetY+pos[0]  , offsetY+pos[0]+S[0]   ), slice(  offsetX+pos[1]  , offsetX+pos[1]+S[1] ) )
        M=numpy.less(0,M)
        mask[  inset   ] =  (num+1)*M
    return mask


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

_cross_data = "\
\x00\x00\x02\x71\
\x89\
\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d\x49\x48\x44\x52\x00\
\x00\x00\x10\x00\x00\x00\x10\x08\x06\x00\x00\x00\x1f\xf3\xff\x61\
\x00\x00\x00\x06\x62\x4b\x47\x44\x00\xff\x00\xff\x00\xff\xa0\xbd\
\xa7\x93\x00\x00\x00\x09\x70\x48\x59\x73\x00\x00\x0d\xd7\x00\x00\
\x0d\xd7\x01\x42\x28\x9b\x78\x00\x00\x00\x07\x74\x49\x4d\x45\x07\
\xde\x07\x10\x0d\x3b\x14\x28\x69\xf4\x66\x00\x00\x01\xfe\x49\x44\
\x41\x54\x38\xcb\x85\x93\xcf\x6a\x14\x41\x10\xc6\x7f\x5d\xdd\x59\
\x77\x3d\x86\x44\x14\xd4\xa0\x41\xfc\x17\xc5\x9b\x20\x1e\x66\x86\
\x3d\xa9\x2f\x10\x04\x4f\x7a\x51\x1f\xc3\x47\x10\x05\x51\x3c\x08\
\x82\x4f\x90\xcb\x30\x33\x37\x11\xd4\x8d\xa2\x17\x49\x02\x1a\x3c\
\x18\x04\x0f\x2e\xe8\x6c\xa6\xbb\x3c\xa4\x27\x4c\x16\xc5\x86\xa6\
\xa1\xaa\xbe\xaf\xbf\xae\xfe\xca\x0c\x87\x49\x8f\x9d\x25\x40\x88\
\xa7\xc6\xfd\xbf\x98\x97\x08\xb6\x31\xe9\x00\x35\x46\x4d\x9e\x57\
\x13\xdf\xc8\x29\xdf\xc8\x25\xef\x65\x2e\xcf\xab\x7a\xba\x0e\xb0\
\xd2\x61\xb4\x80\xcf\xf3\xaa\x6e\xb6\xed\xbd\x34\xc9\x14\x18\x01\
\x25\xca\x66\x9a\x64\x3f\x7c\x23\x57\x22\x91\x6f\x55\x99\xf8\x04\
\x07\x78\x2c\xea\x6b\x59\x03\x8e\xf0\xef\xf5\xb8\xac\x8a\x5b\xc3\
\x61\x32\x00\x42\xab\xc0\xe7\x79\x55\xfb\x5a\x56\x80\x59\xeb\xc2\
\x00\xb8\xb3\x07\x66\xb8\x0c\x9c\x00\x6e\xa6\x49\xb6\x1c\x55\x9b\
\xb6\x39\x92\x26\xd9\x19\x20\x03\xdc\xb1\xa5\xad\xba\xbf\x7f\xfb\
\x51\x4b\x62\x0c\x4b\x22\xe1\x2d\x30\x89\x74\x4f\xe3\x53\x42\x4b\
\x00\x70\x3d\x9e\xbd\xb5\xd5\x83\x1b\x46\x34\x64\xcb\x1f\x1e\x18\
\xd1\x79\x23\x61\x23\x78\x39\x04\x7c\x8e\x35\xfb\xd2\x24\x3b\x07\
\x88\x00\xe2\x7a\x7e\x02\x1c\xdf\x15\x0b\x0b\xbf\xc6\xbd\xd1\xfa\
\xea\x01\x23\xa2\x3f\x83\x97\x59\x60\x7d\xaa\x17\x8b\x80\x0a\x10\
\xfc\xb6\xed\x01\x5b\x9d\xa4\x07\xae\x6d\x7e\x9a\xeb\xe7\x79\x55\
\x97\x55\xf1\x15\x78\x38\x45\x30\x6e\xcd\x23\xaa\x04\xe0\x7d\x4c\
\xd4\xc6\xe8\x61\x6b\xc3\xf7\x10\xe4\x74\x9a\x64\xef\x00\xca\xaa\
\xb8\x0d\x3c\x6b\xd1\xbd\x41\x53\xec\x69\x62\x59\x15\x4f\xe2\xcd\
\xa3\xa2\x2c\xbf\x79\x2f\x17\x50\xde\x00\xe7\xd3\x24\x7b\x1d\x49\
\x6e\x44\xfc\x0b\x0d\xc6\x01\xa1\xf5\x81\x45\x51\x1f\xe4\x22\x4a\
\x05\x7c\x01\x8e\x4e\x49\xfe\x18\x2f\x58\xb4\x2e\xcc\x77\xbf\x71\
\xc7\x89\x06\x71\xd6\xbf\x04\xd2\xbf\x80\x01\xce\x02\x7d\x23\xba\
\x80\xee\xce\x84\xee\xb1\xb1\x62\xac\x75\xe1\x55\x59\x15\x06\xb8\
\x0b\x3c\x07\x72\xe0\x3e\x70\xb5\xac\x8a\x93\x22\x3a\xc6\x60\x76\
\xfd\x35\x1c\x26\xfd\xce\x94\xd1\x25\xb4\x2e\x04\x37\xe3\xc3\xe4\
\xb7\x13\x55\x33\xd3\x99\x81\xb6\xce\xb8\x18\xb4\x31\xd0\xfa\xa0\
\x01\xc4\x37\x22\xbe\x11\xe9\xc6\xa6\xea\xfc\x1f\xdb\x37\xde\x59\
\x25\x68\xce\x04\x00\x00\x00\x00\x49\x45\x4e\x44\xae\x42\x60\x82\
\
"

if __name__=="__main__":
    app=Qt.QApplication([])
    w = mainwindow()
    w.show()
    app.exec_()

