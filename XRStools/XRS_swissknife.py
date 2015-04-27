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

inputtext=""

def main():
    global  inputtext
    filename = sys.argv[1]

    inputfile = open(filename,"r")
    yamlData = load(inputfile, Loader=Loader)
    inputtext = open(filename,"r").read()

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





def load_scans(yamlData):
    """
    This command harvest the selected signals.
    the instructions on the scans to be taken must be in the form( as submembers ofload_scans ) :


    roiaddress :  "hdf5filename:nameofroigroup"  # the same given in create_rois
    expdata    :  "absolutepathtoaspecfile"  # this points to a spec file

    elastic_scans    : [623]
    fine_scans       : [626,630,634,638,642]
    n_loop           : 4
    long_scan        : 624

    signaladdress : "nameofsignalgroup"  # Target group fro writing Relative to ROI (and in the same file)!!!!

    #############################################################
    # OPTIONALS
    #
    order : [0,1,2,3,4,5]  #  list of integers (0-5) which describes the order of modules in which the 
			   #	ROIs were defined (default is VD, VU, VB, HR, HL, HB; i.e. [0,1,2,3,4,5])

    rvd : -41              # mean tth angle of HL module (default is 0.0)
    rvu : 85               # mean tth angle of HR module (default is 0.0)
    rvb : 121.8            # mean tth angle of HB module (default is 0.0)
    rhl : 41.0             # mean tth angle of VD module (default is 0.0)
    rhr : 41.0             # mean tth angle of VU module (default is 0.0)
    rhb : 121.8            # mean tth angle of VB module (default is 0.0)
    #
    """
    mydata = yamlData["load_scans"]
    if mydata is not None and mydata.has_key("active") :
        if mydata["active"]==0:
            return


    roiaddress=None
    roiaddress = mydata["roiaddress"]
    pos = roiaddress.rfind(":")
    if ( pos==-1):
        raise Exception, """
roiaddress   must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""


    filename, groupname = roiaddress[:pos], roiaddress[pos+1:]
    file= h5py.File(filename,"r")
    rois = {}
    shape=xrs_rois.load_rois_fromh5(file[groupname],rois)
    file.close()

    roiob = xrs_rois.roi_object()
    roiob.load_rois_fromMasksDict(rois ,  newshape = shape, kind="zoom")
    
    reader = xrs_read.read_id20(mydata["expdata"] , monitorcolumn='kapraman')


    reader.set_roiObj(roiob)

    elastic_scans   =  mydata["elastic_scans"][:] 
    fine_scans      =  mydata["fine_scans"][:] 
    n_loop          =  mydata["n_loop"] 
    long_scan       =  mydata["long_scan"]

    reader.loadelasticdirect(elastic_scans)
    reader.loadloopdirect(fine_scans,n_loop)
    print " LUNGO " 
    reader.loadlongdirect(long_scan)

    

    reader.getspectrum()

    reader.geteloss()
    reader.gettths(
        rvd   = gvord(mydata,"rvd",0.0) ,
        rvu   = gvord(mydata,"rvu",0.0) ,
        rvb   = gvord(mydata,"rvb",0.0) ,
        rhl   = gvord(mydata,"rhl",0.0) ,
        rhr   = gvord(mydata,"rhr",0.0) ,
        rhb   = gvord(mydata,"rhb",0.0) ,
        order = gvord(mydata,"order",        [0,1,2,3,4,5])
        )

    reader.save_state_hdf5( filename, groupname+"/"+ mydata["signaladdress"], comment = inputtext )

#det value or default
def gvord(yamldata, key, default=None):
    if yamldata.has_key(key):
        return yamldata[key]
    else:
        return default



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

    roiaddress : "myfile.hdf5:/path/to/hdf5/group"

    so that writing will be done in file myfile.hdf5 at the group /path/to/hdf5/group
 
    """
    mydata = yamlData["create_rois"]
    if mydata is not None and mydata.has_key("active") :
        if mydata["active"]==0:
            return
    image4roi = None
    if mydata is not None and mydata.has_key("expdata") :
        repertorio = mydata["expdata"]
        scans = mydata["scans"]
        experiment = xrs_read.read_id20( repertorio ,monitorcolumn='kapraman')
        image4roi =  experiment.SumDirect( scans )

    roiaddress=None
    if mydata is not None and mydata.has_key("roiaddress") :
        roiaddress = mydata["roiaddress"]
        pos = roiaddress.rfind(":")
        if ( pos==-1):
            raise Exception, """
roiaddress for create_rois must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""
        roiaddress = roiaddress[:pos], roiaddress[pos+1:]
    

    from  PyQt4 import Qt, QtCore
    from XRStools import roiSelectionWidget

    app=Qt.QApplication([])
    w4r = roiSelectionWidget.mainwindow()
    if image4roi is not None:
        w4r.showImage( image4roi , xrs_rois.get_geo_informations(image4roi.shape) )
    w4r.show()
    app.exec_()

    if w4r.isOK:
        if roiaddress is  None:
            roiaddress =  Qt.QFileDialog.getSaveFileName()
            roiaddress=str(roiaddress)
            if len(roiaddress)==0: roiaddress=None
            if roiaddress is not None :
                pos = roiaddress.rfind(":")
                if ( pos==-1):
                    raise Exception, """
roiaddress for create_rois must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""
                roiaddress = roiaddress[:pos], roiaddress[pos+1:]
    
            
        if roiaddress is not None:

            h5=h5py.File(roiaddress[0],'a')
            h5.require_group(roiaddress[1]+"/rois_definition")
            masks = w4r.getMasksDict()
            image4roi = w4r.image
            h5group =  h5[roiaddress[1]+"/rois_definition"]
            if(   "image" in h5group.keys() or "roi_dict" in  h5group.keys()):
                raise Exception, (" Rois data already present in  " + roiaddress[0]+":"+roiaddress[1])

            h5group["image"]=image4roi
            h5group.require_group("rois_dict")
            h5group=h5group["rois_dict"]
            xrs_rois.write_rois_toh5(h5group, w4r.getMasksDict() )
            h5.flush()
            h5.close()

             

swissknife_operations={
    "help"        :  help,
    "create_rois" :  create_rois,
    "load_scans" :  load_scans,
}

if __name__=="__main__":
    main()
