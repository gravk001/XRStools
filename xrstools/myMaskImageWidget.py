from  PyMca5.PyMcaGui  import MaskImageWidget  as sole_MaskImageWidget
from  PyQt4 import Qt, QtCore
import numpy 

class MaskImageWidget(sole_MaskImageWidget.MaskImageWidget):
    changeTagOn=False
    def _graphSignal(self, ddict, ownsignal = None):
 
        if self.changeTagOn:
            if self.__selectionMask is not None and ddict['event']=="mouseClicked" and ddict['button']=="middle" :
                x,y = int(ddict["x"]), int(ddict["y"])

                y, x = sole_MaskImageWidget.convertToRowAndColumn(x , y,
                                             self.__imageData.shape,
                                             xScale=self._xScale,
                                             yScale=self._yScale,
                                             safe=True)

                id_target = self.__selectionMask[y,x]
                if id_target and id_target!= self._nRoi:
                    mask_target =  (self.__selectionMask == id_target  )
                    mask_swap   =  (self.__selectionMask== self._nRoi)
                    self.__selectionMask [mask_target] = self._nRoi
                    self.__selectionMask [mask_swap] = id_target
                    emitsignal = True
                    if emitsignal:
                        self.plotImage(update = False)
                        self._emitMaskChangedSignal()
                    return
        # super(MaskImageWidget, self)._graphSignal(ddict, ownsignal) 
        sole_MaskImageWidget.MaskImageWidget._graphSignal(self, ddict, ownsignal)


    def dragEnterEvent(self,event):
        print dir(event.mimeData())
        print list(event.mimeData().formats())
        print event.mimeData().text()

        model = Qt.QStandardItemModel()
        model.dropMimeData(event.mimeData(), QtCore.Qt.CopyAction, 0,0, Qt.QModelIndex())
        print model

        if event.mimeData().hasFormat('application/x-qabstractitemmodeldatalist'):
            event.acceptProposedAction()

            bytearray = event.mimeData().data('application/x-qabstractitemmodeldatalist')
            data_items = decode_data(bytearray)
            print data_items

        if event.mimeData().hasFormat("text/plain"):
            print " OK "
            event.acceptProposedAction()
 
    def dropEvent(self, e):

        position = e.pos()    
        print "POSITION ", position.x(), position.y()
        x,y = position.x(), position.y()
        
        x,y = self.graph.displayXY2dataXY(x,y ) 
    
        mask = self.getSelectionMask()
        
        Ct = mask[int(y), int(x)]
        if Ct:
            bytearray = e.mimeData().data('application/x-qabstractitemmodeldatalist')
            data = decode_data(bytearray)
            Cc = (2-data[0])*4+data[1]
            zonet = (mask==Ct)
            zonec = (mask==Cc)
            mask[zonet]=Cc
            mask[zonec]=Ct

            self.annotateSpots()


    def annotateSpots(self):

        self.graph.clearMarkers()
        mask = self.getSelectionMask().astype("i")

        nspots = mask.max()
        for i in range(1,nspots+1):
            m = (mask==i).astype("f")
            msum=m.sum()
            if msum:
                ny,nx = m.shape
                px= (m.sum(axis=0)*numpy.arange(nx)).sum()/msum
                py= (m.sum(axis=1)*numpy.arange(ny)).sum()/msum

                print " ##################################    ", px, py

                res=self.graph.insertMarker( px, py, "legend", "%d"%(i), color='black', selectable=False, draggable=False, searchFeature=True, 
                                           xytext = (-20, 0),
                                           textcoords = 'offset points',
                                           ha = 'right', va = 'bottom',
                                           bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.4),
                                           arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))




def decode_data( bytearray):
        
    data = []
    item = {}
    
    ds = QtCore.QDataStream(bytearray)
    while not ds.atEnd():
            
        row = ds.readInt32()
        column = ds.readInt32()

        return row, column

        print row, column

        map_items = ds.readInt32()
        for i in range(map_items):
               
            key = ds.readInt32()
            
            value = QtCore.QVariant()
            ds >> value
            item[QtCore.Qt.ItemDataRole(key)] = value
               
        data.append(item)
           
    return data
