import string
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

    the instructions on the scans to be taken must be in the form( as submembers of create_rois) :

    expdata :  "absolutepathtoaspecfile"  # this points to a spec file
    scans   : [623,624]                   # a list containing one or more scans 
                                          # for the elastic part.
                                          # They will be summed to an image
                                          # and rois will then be extracted from this image.
    
    If expdata is not given in the yaml input file, then the roicreation widget will be launched
    without data, and you have to open some file for the image.

    Before quitting the roiwidget you have the possibility to write the rois into a hdf5 file.
    You do this by choosing a filename for the hdf5 file, and a name for the roi selection.
    An entry, with this matter name, will be created into the hdf5 file ( all subsequent
    treatements, done by other methods than this, which uses this roi selection 
    ( called by its name )  will be reported into subentries of such entry)

    You can also decide the target entry already in the yaml file, using such from :

    target : "myfile.hdf5:/path/to/hdf5/group"

    so that writing will be done in file myfile.hdf5 at the group /path/to/hdf5/group
 
    """
    mydata = yamlData["create_rois"]
    image4roi = None
    if mydata is not None and mydata.has_key("expdata") :
        repertorio = mydata["expdata"]
        scans = mydata["scans"]
        experiment = xrs_read.read_id20( repertorio ,monitorcolumn='kapraman')
        image4roi =  experiment.SumDirect( scans )

    target=None
    if mydata is not None and mydata.has_key("target") :
        target = mydata["target"]
        pos = target.rfind(":")
        if ( pos==-1):
            raise Exception, """
target for create_rois must be given in the form  target : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""
        target = target[:pos], target[pos+1:]
    

    from  PyQt4 import Qt, QtCore
    from XRStools import roiSelectionWidget

    app=Qt.QApplication([])
    w4r = roiSelectionWidget.mainwindow()
    if image4roi is not None:
        w4r.showImage( image4roi , xrs_rois.get_geo_informations(image4roi.shape) )
    w4r.show()
    app.exec_()

    if w4r.isOK:
        if target is  None:
            target =  Qt.QFileDialog.getSaveFileName()
            target=str(target)
            if len(target)==0: target=None
            if target is not None :
                pos = target.rfind(":")
                if ( pos==-1):
                    raise Exception, """
target for create_rois must be given in the form  target : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""
                target = target[:pos], target[pos+1:]
    
            
            

        if target is not None:

            h5=h5py.File(target[0],'a')
            h5.require_group(target[1]+"/rois_definition")
            masks = w4r.getMasksDict()
            image4roi = w4r.image
            h5group =  h5[target[1]+"/rois_definition"]
            if(   "image" in h5group.keys() or "roi_dict" in  h5group.keys()):
                raise Exception, (" Rois data already present in  " + target[0]+":"+target[1])

            h5group["image"]=image4roi
            h5group.require_group("rois_dict")
            h5group=h5group["rois_dict"]
            write_rois_toh5(h5group, w4r.getMasksDict() )
            h5.flush()
            h5.close()

             

def  write_rois_toh5(h5group,md):
    for key in md.keys():
        h5group.require_group(key)
        h5group[key]["origin"]=md[key][0]
        h5group[key]["mask"]=md[key][1]
    

swissknife_operations={
    "help"        :  help,
    "create_rois" :  create_rois,
}

if __name__=="__main__":
    main()
