__doc__="""
Documentation of XRS_swissknife
-------------------------------

The script is called in this way ::

  XRS_swissknife yourinput.yaml

The input file is in *yaml* format. In *yaml* format each line introduces an item
and the indentation expresses the hierarchy.
An example is ::

  Fancy_Reduction :
        parameterTom : 3.14
        parameterJerry : False
        choicesBob     : [1,2,3]

In this example we create an item called *Fancy_Reduction* which contains three subitems.

The *XRS_swissknife* expects that for each operation that you want to apply, you provide an item
named as the operation, and the associated subitems that provide that values for the parameters.

*XRS_swissknife* will do what you want provided you respect the proper indentation. A thing which helps
is using emacs and activating the *python mode*, because python uses the same indentation principle
to structure the code. 

Each processing item has the additional, optional, key **active**.
This key can be set to **0** or **1** to desactivate or not (default is **1**, active) the processing.
Here a desactivation example ::

  Fancy_Reduction :
        active : 0
        parameterTom : 3.14
        parameterJerry : False
        choicesBob     : [1,2,3]

The following documentation has been created automatically, for each functionality, based on the documentation string
written in the code for the fonctionality itself.  You can also write mathematical expressions :  :math:`\\int{ x dx} = \\frac { x^2 }{ 2 }`
and even include graphics.

"""
import string
import numpy as np
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper            

import h5py

import sys
import xrs_rois
import xrs_read
import theory
import extraction
import xrs_prediction
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
    **help**

    Displays doc on the operations. In the input file ::

       help :

    will trigger printing of all the available operation names ::

       help :
           create_rois
           load_scans

    will print the help on *create_rois* and the help about *load_scans*.
    By the way, it is the same that you can read here because the *sphinx* doc-generation tool
    reads the same docstrings contained in the code.
     
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


def Extraction(yamlData):
    """ 
    **Extraction**

    Function to extract the interesting signal after removal of Compton profile,
    linear baselines,  Pearson profile
    Example ::

     Extraction :
         active : 1
         dataadress   : "pippo.h5:/ROI_A/loaded_datas"         # where load_scans wrote data
         HFaddress    : "pippo.h5:/ROI_A/loaded_datas/HF_O"    # where compton profiles have been calculated
         # prenormrange : [ 5 , .inf ]		

         analyzerAverage :                                     #  averaging over analysers
             active : 1
             which : [0,11  , 36,59   ]
             errorweighing  : False

         removeLinearAv :                                      #  fit a linear baseline and remove it
             active  : 1
             region1 :  [520.0,532.0]   
             region2 :  None 
             ewindow : 100 
             scale : 1

         removePearsonAv:                                      # fit a Pearson and remove it
             active  : 0
             region1 :  [520.0,532.0]   
             region2 :  None  
             guess :
                 Peak_position : 600.0
                 FWHM          :  10
                 Shape         : "Lorentzian" 
                 Peak_intensity: 100.0
                 linear_slope  : 1
                 linear_background : 0
                 scaling_factor : 1

         view   :   0
         target :   "myextraction"                            # path relative to dataadress where extracted signal will be written

    """

    mydata = yamlData["Extraction"]
    if  mydata.has_key("active") :
        if mydata["active"]==0:
            return   
    reader , filename, groupname= read_reader(mydata, name="dataadress")
    HF     = read_HF(mydata, name="hfspectrum_address")

    extr  = extraction.extraction(reader , HF)

    if mydata.has_key("analyzerAverage"):
        aa_data = mydata["analyzerAverage"]
        if gvord(aa_data,"active",True):
            which = aa_data["which"]
            errorweighing  = gvord(aa_data,"errorweighing",False)
            extr .analyzerAverage(which,errorweighing=errorweighing)

    if mydata.has_key("removeLinearAv"):
        rla_data = mydata["removeLinearAv"]
        if gvord(rla_data,"active",True):
            region1 = rla_data["region1"]
            region2 = gvord( rla_data,"region2",None)
            ewindow = gvord( rla_data,"ewindow",100)
            scale = gvord( rla_data,"scale",100)
            extr .removeLinearAv(region1, region2=region2,ewindow=ewindow, 
                                 scale=scale, view = gvord(mydata,"view",False),                                 
                             ) 
    print gvord(mydata,"view",False)
  
    extr.save_state_hdf5( filename, groupname+"/"+ mydata["target"]+"/"+"datas", comment = inputtext )
   


def read_HF(mydata, name="hfspectrum_address"):

    dataadress = mydata[name]
    pos = dataadress.rfind(":")
    if ( pos==-1):
        raise Exception, """
hfspectrum_address   must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""
    filename, groupname = dataadress[:pos],dataadress [pos+1:]
    HF = theory.HFspectrum(None,None,None, initialise=False)
    HF.load_state_hdf5( filename, groupname+"/"+"datas")
    return HF


def  HFspectrum(yamlData):
    """
    **HFspectrum**

    function for building S(q,w) from tabulated Hartree-Fock Compton profiles to use in
    the extraction algorithm.

    EXAMPLE ::

      dataadress : "hdf5filename:full_nameof_signals_group"   # where load_scans wrote data
      formulas   :  ['O']     # list of strings of chemical sum formulas of which the sample is made up
      concentrations : [1.0]  # list of concentrations of how the different chemical formulas are mixed (sum should be 1)
      correctasym    : [[0.0,0.0,0.0]]  #  single value or list of scaling values for the HR-correction to 
          		              # the 1s, 2s, and 2p shells. one value per element in the list of formulas
      hfspectrum_address : "nameofgroup" # Target group for writing Relative to dataadress (and in the same file)!!!!

 
    """

    mydata = yamlData["HFspectrum"]
    if mydata is not None and mydata.has_key("active") :
        if mydata["active"]==0:
            return   

    reader , filename, groupname= read_reader(mydata, name="dataadress")

    hf   = theory.HFspectrum(reader ,   
                             mydata["formulas"]     ,
                             mydata["concentrations"]     ,
                             mydata["correctasym"]  
                             )
    hf.save_state_hdf5( filename, groupname+"/"+ mydata["hfspectrum_address"]+"/"+"datas", comment = inputtext )
    

def load_scans(yamlData):
    """
    **load_scans**

    This command harvest the selected signals.
    the instructions on the scans to be taken must be in the form( as submembers ofload_scans ) ::


     load_scans :
         roiaddress :  "hdf5filename:nameofroigroup"  # the same given in create_rois
         expdata    :  "absolutepathtoaspecfile"  # this points to a spec file

         elastic_scans    : [623]
         fine_scans       : [626,630,634,638,642]
         n_loop           : 4
         long_scan        : 624

         signaladdress : "nameofsignalgroup"  # Target group for writing Relative to ROI (and in the same file)!!!!

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

    reader.save_state_hdf5( filename, groupname+"/"+ mydata["signaladdress"]+"/"+"datas", comment = inputtext )

#det value or default
def gvord(yamldata, key, default=None):
    if yamldata.has_key(key):
        return yamldata[key]
    else:
        return default



def create_rois(yamlData):
    """
    **create_rois**

    This launches the roi selector.
    If yamlData contains instructions on the scan to be taken , the roi 
    selector will pop up with the image ready to be used. Otherwise you have to
    browse some file around.

    At the end you have the possibility to write the rois on a container hdf5 file.

    In the extreme case when you give no argument ( parameters) ::

        create_rois :

    The roi selector window  pops up and you have to select an image at the beginning
    and select a hdf5 file, and a name of a node inside the file, where the rois 
    will be written.

    example ::

      create_rois :

        expdata :  "absolutepathtoaspecfile"  # this points to a spec file
        scans   : [623,624]                   # a list containing one or more scans 
                                              # for the elastic part.
                                              # They will be summed to an image
                                              # and rois will then be extracted from this image.
        roiaddress : "myfile.hdf5:/path/to/hdf5/group"  # the target destination for rois
   
    If expdata is not given in the yaml input file, then the roicreation widget will be launched
    without data, and you have to open some file for the image.

    If roiaddress is not given,  before quitting the roiwidget you have the possibility to write the rois into a hdf5 file.
    You do this by choosing a filename for the hdf5 file, and a name for the roi selection.
    An entry, with this matter name, will be created into the hdf5 file ( all subsequent
    treatements, done by other methods than this, which uses this roi selection 
    ( called by its name )  will be reported into subentries of such entry)

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


def XRSprediction(yamlData):
	""" **prediction**

	This launches the XRS prediction routines.

	If yamlData contains information about: the sample, the incident beam, 
	the analyzer, the detector, the polarization, and the HF compton profiles,
	this will create the desired predicted XRS data.

	At the end you have the possibility to write the predicted profiles into a container hdf5 file.

	In the extreme case when you give no argument ( parameters) ::

		xrs_prediction :

	The following canonical example will be run.

	example ::

		xrs_prediction :
			active : 1

			sample :
				chem_formulas : ['C']	# list of strings of chemical sum formulas
				concentrations : [1.0]	# list of concentrations, should contain values between 0.0 and 1.0
				densities : [2.266]		# list of densities of the constituents [g/cm^3]
				angle_tth : 35.0		# scattering angle [deg]
				sample_thickness : 0.1	# sample thickness/diameter in [cm]
				angle_in : None			# incident beam angle in [deg] relative to sample surface normal
				angle_out : None		# beam exit angle in [deg] relatice to sample surface normal
										# (negative for transmission geometry)
				shape : 'sphere'		# keyword, can be 'slab' or 'sphere'
				molar_masses : [12.0]	# list of molar masses of all constituents

			incident_beam :
				i0_intensity : 1e13		#  # number of incident photons [1/sec]
				beam_height : 10.0		# in micron
				beam_width : 20.0		# in micron
				
			analyzer :
				material : 'Si'			# analyzer material (e.g. 'Si', 'Ge')
				hkl : [6,6,0]			# [hkl] indices of reflection used
				mask_d : 60.0			# analyzer mask diameter in [mm]
				bend_r : 1.0			# bending radius of the crystal [mm]
				energy_resolution : 0.5 # energy resolution [eV]
				diced : False			# boolean (True or False) if a diced crystal is used or not (defalt is False)
				thickness : 500.0		# thickness of the analyzer crystal
				database_dir : installation_dir

			compton_profiles :
				eloss_range : np.arange(0.0,1000.0,0.1)
				E0 : 9.7

			detector :
				energy : 9.7			# analyzer energy [keV]
				thickness : 500.0		# thickness of the active material [microns]
				material : 'Si'			# detector active material

			thomson :
				scattering_plane : 'vertical'	# keyword to indicate scattering plane relative to lab frame ('vertical' or 'horizontal')
				polarization : 0.99				# degree of polarization (close to 1.0 for undulator radiation)

			saveaddress : "myfile.hdf5:/path/to/hdf5/group"  # the target destination, if data should be saved
   
	"""
	mydata = yamlData["XRSprediction"]
	if mydata is not None and mydata.has_key("active"):
		if mydata["active"]==0:
			return

	if mydata is not None:
		try:
			sample = mydata["sample"]
		except:
			sample = {}
		chem_formulas    = gvord(sample,"chem_formulas",["C"])
		concentrations   = gvord(sample,"concentrations", [1.0])
		densities        = gvord(sample,"densities", [2.266])
		angle_tth        = gvord(sample,"angle_tth", 35.0)
		sample_thickness = gvord(sample,"sample_thickness", 0.1)
		angle_in         = gvord(sample,"angle_in", None)
		angle_out        = gvord(sample,"angle_out", None)
		shape            = gvord(sample,"shape", 'sphere')
		molar_masses     = gvord(sample,"molar_masses", [12.0])
		sample_obj = xrs_prediction.sample(chem_formulas, concentrations, densities, angle_tth, sample_thickness, angle_in, angle_out, shape, molar_masses)
	
		try:
			incident_beam = mydata["incident_beam"]
		except:
			incident_beam = {}
		i0_intensity  = gvord(incident_beam,"i0_intensity", 1e13)
		beam_height   = gvord(incident_beam,"beam_height", 10.0)
		beam_width    = gvord(incident_beam,"beam_width", 20.0)
		beam_obj = xrs_prediction.beam(i0_intensity, beam_height, beam_width, 0.0)

		try:
			analyzer = mydata["analyzer"]
		except:
			analyzer = {}
		material          = gvord(analyzer,"material", 'Si')
		hkl               = gvord(analyzer,"hkl", [6,6,0])
		mask_d            = gvord(analyzer,"mask_d", 60.0)
		bend_r            = gvord(analyzer,"bend_r", 1.0)
		energy_resolution = gvord(analyzer,"energy_resolution", 0.5)
		diced             = gvord(analyzer,"diced", False)
		thickness         = gvord(analyzer,"thickness", 500.0)
		database_dir      = gvord(analyzer,"database_dir", "/home/christoph/sources/XRStools/data/chitables/")
		analyzer_obj = xrs_prediction.analyzer(material, hkl, mask_d, bend_r, energy_resolution, diced, thickness, database_dir)

		try:
			detector = mydata["detector"]
		except:
			detector = {}
		energy    = gvord(detector,"energy", 9.7)
		thickness = gvord(detector,"thickness", 500.0)
		material  = gvord(detector,"material", 'Si')
		detector_obj  = xrs_prediction.detector(energy, thickness, material, [256,768])

		try:
			compton_profile = mydata["compton_profile"]
		except:
			compton_profile = {}
		eloss_range = gvord(compton_profile,"eloss_range", np.arange(0.0,1000.0,0.1))
		E0          = gvord(compton_profile,"E0", 9.7)
		compton_profile_obj = xrs_prediction.compton_profiles(sample_obj, eloss_range, E0)

		try:
			thomson = mydata["polarization"]
		except:
			thomson = {}
		scattering_plane = gvord(thomson,"scattering_plane", 'vertical')
		polarization     = gvord(thomson,"polarization", 0.99)
		thomson_obj = xrs_prediction.thomson(compton_profile_obj.get_energy_in_keV(),compton_profile_obj.get_E0(),compton_profile_obj.get_tth())

		abs_cross_section_obj = xrs_prediction.absolute_cross_section(beam_obj, sample_obj, analyzer_obj, detector_obj, thomson_obj, compton_profile_obj)
		abs_cross_section_obj.plot_abs_cross_section()


def read_reader(mydata, name="dataadress"):

    dataadress = mydata[name]
    pos = dataadress.rfind(":")
    if ( pos==-1):
        raise Exception, """
dataaddress   must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
"""
    filename, groupname = dataadress[:pos],dataadress [pos+1:]
    reader = xrs_read.read_id20(None)
    reader.load_state_hdf5( filename, groupname+"/"+"datas")
    return reader, filename, groupname

  

swissknife_operations={
	"help"        :  help,
	"create_rois" :  create_rois,
	"load_scans" :  load_scans,
	"HFspectrum"  : HFspectrum,
	"Extraction"    :  Extraction,
	"XRSprediction" : XRSprediction
}

if __name__=="__main__":
    main()

