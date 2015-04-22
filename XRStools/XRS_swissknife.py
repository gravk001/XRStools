from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper            

import h5py

import sys
import xrs_rois
import xrs_read


def main():
    filename = sys.argv[1]

    yamlData = load(open(filename,"r"), Loader=Loader)

    for key in yamlData.keys():
        swissknife_operations[key](yamlData   )
    
    # workbench_file = yamlData["workbench_file"]


def help(yamlData):
    """
    Displays all operations and their doc
    """
    print " HELP " *15
    if yamlData ["help"] is None:
        print """
              Printing all the function names
              To get help on a specific function:

              help : "functionName"
"""
        for key,func in swissknife_operations.iteritems():
            print " FUNCTION : " , key

    else:
        func = swissknife_operations[ yamlData ["help"]]
        print "---------------------------------------"
        print func.__doc__

def create_rois(yamlData):
    """
    This launches the roi selector.
    If yamlData contains instructions on the scan to be taken , the roi 
    selector will pop up with the image ready to be used. Otherwise you have to
    browse some file around.

    At the end you have the possibility to write the rois on a container hdf5 file

    the instructions on the scans to be taken must be in the form :

    expdata :  "absolutepathtoaspecfile"  # this points to a spec file
    scans   : [623,624]                   # a list containing one or more scans 
                                          # for the elastic part.
                                          # They will be summed to an image
                                          # and rois will then be extracted from this image.
    
    If expdata is not given in the yaml input file, then the roicreation widget will be launched
    without data, and you have to open some file for the image.

    Before quitting the roiwidget you have the possibility to write the rois into a hdf5 file.
    You do this by choosing a filename for the hdf5 file, and a name for the roi selection.
    An entry, with this matter name, will be created into the hdf5 file, and all subsequent
    treatements which uses this roi selection ( called by its name )  will be reported into subentries
 
    """
    mydata = yamlData["create_rois"]
    image4roi = None
    if mydata is not None and mydata.has_key("expdata") :
        repertorio = mydata["expdata"]
        scans = mydata["scans"]
        experiment = xrs_read.read_id20( repertorio ,monitorcolumn='kapraman')
        image4roi =  experiment.SumDirect( scans )

    from  PyQt4 import Qt, QtCore
    from XRStools import roiSelectionWidget

    app=Qt.QApplication([])
    w4r = roiSelectionWidget.mainwindow()
    if image4roi is not None:
        w4r.showImage( image4roi , xrs_rois.get_geo_informations(image4roi.shape) )
    w4r.show()
    app.exec_()

    


swissknife_operations={
    "help"        :  help,
    "create_rois" :  create_rois,
}

if __name__=="__main__":
    main()
