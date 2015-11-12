import numpy as np
import matplotlib.pyplot as plt

from matplotlib.path import Path
import  xrs_utilities 
import math_functions 
from matplotlib.widgets import Cursor, Button
from scipy.ndimage import measurements
from scipy import signal
import  xrs_rois 
import xrs_scans
		

def findroisColumns(scans,scannumbers,roi_obj, whichroi,logscaling=False):
        """
        Constructs a waterfall plot from the columns in a scan, i.e. energy vs. pixel number along the ROI
        scannumbers = scannumber or list of scannumbers from which to construct the plot
        whichroi    = integer (starting from 0) from which ROI to use for constructing the plot
        """
        if not isinstance(scannumbers,list):
                scannums = []
                scannums.append(scannumbers)
        else:
                scannums = scannumbers 

        if not roi_obj.indices:
                'Please define some zoom ROIs first.'
                return
        if not roi_obj.kind == 'zoom':
                'Currently this feature only works for ROIs of type \'zoom\'.'
                return

        xinds     = np.unique(roi_obj.x_indices[whichroi])
        yinds     = np.unique(roi_obj.y_indices[whichroi])
        scanname  = 'Scan%03d' % scannums[0]
        edfmats   = np.zeros_like(scans[scanname].edfmats)
        energy    = scans[scanname].energy
        waterfall = np.zeros((len(yinds),len(energy)))

        for scannum in scannums:
                scanname = 'Scan%03d' % scannum
                scanedf  = scans[scanname].edfmats
                scanmonitor = scans[scanname].monitor
                for ii in range(len(energy)):
                        edfmats[ii,:,:] += scanedf[ii,:,:]/scanmonitor[ii]

        for ii in range(len(energy)):
                for jj in range(len(yinds)):
                        waterfall[jj,ii] = np.sum(edfmats[ii,xinds,yinds[jj]])

        plt.figure()
        for ii in range(len(yinds)):
                plt.plot(waterfall[ii,:])

        fig = plt.figure()
        ax = fig.add_subplot(111)

        if logscaling:
                ax.imshow(np.log(np.transpose(waterfall)), interpolation='nearest')
        else:
                ax.imshow(np.transpose(waterfall), interpolation='nearest')

        ax.set_aspect('auto')
        plt.xlabel('ROI pixel')
        plt.ylabel('energy point')
        plt.show()


def get_auto_rois_eachdet(scans,DET_PIXEL_NUM ,scannumbers,kernel_size=5,threshold=100.0,logscaling=True,colormap='jet',interpolation='bilinear'):
    """
    Define ROIs automatically using median filtering and a variable threshold for each detector
    separately.
    scannumbers   = either single scannumber or list of scannumbers
    kernel_size   = used kernel size for the median filter (must be an odd integer)
    logscaling    = set to 'True' if images is to be shown on log-scale (default is True)
    colormap      = string to define the colormap which is to be used for display (anything 
                    supported by matplotlib, 'jet' by default)
    interpolation = interpolation scheme to be used for displaying the image (anything
                    supported by matplotlib, 'nearest' by default)
    """
    # check that kernel_size is odd
    if not kernel_size % 2 == 1:
            print 'The \'kernal_size\' must be an odd number.'
            return

    # create a big image
    image = xrs_scans.create_sum_image(scans,scannumbers)

    # break down the image into 256x256 pixel images
    det_images, offsets = xrs_rois.break_down_det_image(image,DET_PIXEL_NUM)

    # create one roi_object per sub-image
    temp_objs = []
    for ii in range(det_images.shape[0]):
            temp = roi_finder()
            temp.get_auto_rois(det_images[ii,:,:],kernel_size=kernel_size,threshold=threshold,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
            temp_objs.append(temp)

    # merge all roi_objects into one
    merged_obj   = xrs_rois.merge_roi_objects_by_matrix(temp_objs,image.shape,offsets,DET_PIXEL_NUM)
    roi_obj = merged_obj
    return roi_obj

def get_polygon_rois_eachdet(scans,DET_PIXEL_NUM,  scannumbers,logscaling=True,colormap='jet',interpolation='nearest'):
        """
        Define a polygon shaped ROI from an image constructed from 
        the sum of all edf-files in 'scannumbers'
        image_shape = tuple with shape of the current image (i.e. (256,256))
        scannumbers = either single scannumber or list of scannumbers
        logscaling  = set to 'True' if images is to be shown on log-scale (default is True)
        colormap    = string to define the colormap which is to be used for display (anything 
                      supported by matplotlib, 'jet' by default)
        interpolation = interpolation scheme to be used for displaying the image (anything
                        supported by matplotlib, 'nearest' by default)
        """

        # create a big image
        image = xrs_scans.create_sum_image(scans,scannumbers)

        # break down the image into 256x256 pixel images
        det_images, offsets = xrs_rois.break_down_det_image(image,DET_PIXEL_NUM)

        # create one roi_object per sub-image
        temp_objs = []
        for modind in range(det_images.shape[0]):
                temp = roi_finder()
                temp.get_polygon_rois( det_images[modind,:,:],modind,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
                temp_objs.append(temp)

        # merge all roi_objects into one
        merged_obj   = xrs_rois.merge_roi_objects_by_matrix(temp_objs,image.shape,offsets,DET_PIXEL_NUM)
        roi_obj = merged_obj
        return roi_obj

def get_zoom_rois(scans,scannumbers,logscaling=True,colormap='jet',interpolation='nearest'):
        # create a big image
        image = xrs_scans.create_sum_image(scans,scannumbers)

        # create one roi_object per sub-image
        roi_obj = roi_finder()
        roi_obj.get_zoom_rois(image,logscaling=logscaling,colormap=colormap,interpolation=interpolation)        

        return roi_obj.roi_obj



class roi_finder:
	"""
	Class to define ROIs from a 2D image.
	"""

	def __init__(self):
		self.roi_obj = xrs_rois.roi_object() # empty roi object

	def deleterois(self):
		"""
		Clear the existing ROIs by creating a fresh roi_object.
		"""
		self.roi_obj = xrs_rois.roi_object()

	def get_linear_rois(self,input_image,logscaling=True,height=5,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by clicking two points on a 2D image.
		number_of_rois = integer defining how many ROIs should be determined
		input_object   = 2D array, scan_object, or dictionary of scans to define the ROIs from
		logscaling     = boolean, to determine wether the image is shown on a log-scale (default = True)
		height         = integer defining the height (in pixels) of the ROIs
		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# check that the input is a 2d matrix
		if not len(input_image.shape) == 2:
			print 'Please provide a 2D numpy array as input!'
			return

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# save input image for later use
		self.roi_obj.input_image = input_image

		# prepare a figure
		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		# Initialize suptitle, which will be updated
		titlestring = 'Start by clicking the \'Next\'-button.'
		titleInst=plt.suptitle(titlestring)

		# generate an image to be displayed
		figure_obj = plt.imshow(input_image,interpolation=interpolation)

		# set the colormap for the image
		figure_obj.set_cmap(colormap)

		rois   = []
		class Index:
			ind  = 0
			def next(self, event):
				self.ind += 1
				titlestring = 'Click two points for ROI Nr. %02d, hit \'Next\' to continue, \'Finish\' to end.' % self.ind
				titleInst.set_text(titlestring) # Update title
				one_roi = define_lin_roi(height,input_image.shape)
				for index in one_roi:
					input_image[index[0],index[1]] += 1.0e6
				figure_obj.set_data(input_image)
				plt.hold(True)
				rois.append(one_roi)
			def prev(self, event):
				self.ind -= 1
				try:
					titlestring = 'Click the \'Next\' button again to continue.'
					titleInst.set_text(titlestring) # Update title
					# print titlestring
					for index in rois[-1]:
						input_image[index[0],index[1]] -= 1.0e6
					figure_obj.set_data(input_image)
					plt.hold(True)
					rois.pop()
				except:
					pass
			def close(self, event):
				plt.close()
			def dmy(self, event):
				pass # adding a dummy function for the dummy button

		callback = Index()
		axprev   = plt.axes([0.5, 0.05, 0.1, 0.075])
		axnext   = plt.axes([0.61, 0.05, 0.1, 0.075])
		axclose  = plt.axes([0.72, 0.05, 0.1, 0.075])
		axdmy    = plt.axes([0.001, 0.001, 0.001, 0.001]) # for some reason the first botton disappears when clicked
		bdmy     = Button(axdmy,'')                        # which is why I am including a dummy button here
		bdmy.on_clicked(callback.dmy)                   # this way, the real buttons work
		bnext    = Button(axnext, 'Next')
		bnext.on_clicked(callback.next)
		bprev    = Button(axprev, 'Back')
		bprev.on_clicked(callback.prev)
		bclose   = Button(axclose, 'Finish')
		bclose.on_clicked(callback.close)
		plt.show()

		# assign the defined rois to the roi_object class
		self.roi_obj.roi_matrix     = xrs_rois.convert_inds_to_matrix(rois,input_image.shape)
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = rois 
		self.roi_obj.kind           = 'linear'
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(rois)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(rois)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def get_zoom_rois(self,input_image,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by clicking two points on a 2D image.
		number_of_rois = integer defining how many ROIs should be determined
		input_object   = 2D array, scan_object, or dictionary of scans to define the ROIs from
		logscaling     = boolean, to determine wether the image is shown on a log-scale (default = True)
		height         = integer defining the height (in pixels) of the ROIs
		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# check that the input is a 2d matrix
		if not len(input_image.shape) == 2:
			print 'please provide a 2D numpy array as input!'
			return                

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# save input image for later use
		self.roi_obj.input_image = input_image

		# prepare a figure
		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		# Initialize suptitle, which will be updated
		titlestring = 'Start by clicking the \'Next\' button.'
		titleInst=plt.suptitle(titlestring)

		# generate an image to be displayed
		figure_obj = plt.imshow(input_image,interpolation=interpolation)

		# activate the zoom function already
		thismanager = plt.get_current_fig_manager()
		thismanager.toolbar.zoom()

		# set the colormap for the image
		figure_obj.set_cmap(colormap)
		#plt.colorbar()

		# initialize a matrix for the rois (will be filled with areas of ones, twos, etc
		rois = []

		# print info to start: 
		print 'Start by clicking the \'Next\' button.'

		class Index:
			ind          = 0
			initstage    = True
			next_clicked = False
			def next(self, event):
				# for some reason, the first time this is used, it doesn't work, so here is one dummy round
				if self.initstage:
					#self.ind += 1
					self.initstage = False
					plt.sca(ax)
					one_roi = define_zoom_roi(input_image,verbose=True)
					#for index in one_roi:
					#	input_image[index[0],index[1]] *= 1.0
					# reset the matrix to be displayed
					figure_obj.set_data(input_image)
					# reset the zoom
					plt.xlim(0.0,input_image.shape[1])
					plt.ylim(input_image.shape[0],0.0)
					plt.draw()
					titlestring = 'Zoom in to define ROI Nr. %02d, hit \'Next\' to continue.' % (self.ind + 1)
					titleInst.set_text(titlestring) # Update title
				else:
					self.ind += 1
					plt.sca(ax)
					one_roi = define_zoom_roi(input_image)
					for index in one_roi:
						input_image[index[0],index[1]] += 1.0e10
					# reset the matrix to be displayed
					figure_obj.set_data(input_image)
					# reset the zoom
					plt.xlim(0.0,input_image.shape[1])
					plt.ylim(input_image.shape[0],0.0)
					plt.draw()
					rois.append(one_roi)
					titlestring = 'Zoom in to define ROI Nr. %02d, hit \'Next\' to continue, \'Finish\' to end.' % (self.ind + 1)
					titleInst.set_text(titlestring) # Update title
                                                
			def prev(self, event):
				self.ind -= 1
				titlestring = 'Undoing ROI Nr. %02d. Zoom again, click the \'Next\' button to continue.' % (self.ind + 1)
				titleInst.set_text(titlestring) # Update title
				#thedata[roimatrix == self.ind+1]   -= 1.0e6
				#roi_matrix[roimatrix == self.ind+1] = 0.0
				for index in rois[-1]:
					input_image[index[0],index[1]] -= 1.0e10
				figure_obj.set_data(input_image)
				plt.hold(True)
				rois.pop()

			def close(self, event):
				plt.sca(ax)
				one_roi = define_zoom_roi(input_image)
				for index in one_roi:
					input_image[index[0],index[1]] += 1.0e10
				# reset the matrix to be displayed
				figure_obj.set_data(input_image)
				# reset the zoom
				plt.xlim(0.0,input_image.shape[1])
				plt.ylim(input_image.shape[0],0.0)
				plt.draw()
				rois.append(one_roi)
				titlestring = 'Last ROI is Nr. %02d.' % (self.ind + 1)
				titleInst.set_text(titlestring) # Update title
				plt.close()

			def dmy(self, event):
				pass # adding a dummy function for the dummy button

		callback = Index()
		axprev   = plt.axes([0.5, 0.05, 0.1, 0.075])
		axnext   = plt.axes([0.61, 0.05, 0.1, 0.075])
		axclose  = plt.axes([0.72, 0.05, 0.1, 0.075])
		axdmy    = plt.axes([0.001, 0.001, 0.001, 0.001]) # for some reason the first botton disappears when clicked
		bdmy     = Button(axdmy,'')                       # which is why I am including a dummy button here
		bdmy.on_clicked(callback.dmy)                     # this way, the real buttons work
		bnext    = Button(axnext, 'Next')
		bnext.on_clicked(callback.next)
		bprev    = Button(axprev, 'Back')
		bprev.on_clicked(callback.prev)
		bclose   = Button(axclose, 'Finish')
		bclose.on_clicked(callback.close)
		plt.show()

		# assign the defined rois to the roi_object class
		self.roi_obj.roi_matrix     = (xrs_rois.convert_inds_to_matrix(rois,input_image.shape))
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = rois 
		self.roi_obj.kind           = 'zoom'
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(rois)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(rois)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def get_auto_rois(self,input_image,kernel_size=5,threshold=100.0,logscaling=True,colormap='jet',interpolation='bilinear'):
		"""
		Define ROIs by choosing a threshold using a slider bar under the figure. In this function, the entire 
		detector is shown.
		input_image = 2D numpy array with the image to be displayed
		kernal_size = integer defining the median filter window (has to be odd)
		theshold    = initial number defining the upper end value for the slider bar (amax(input_image)/threshold defines this number), can be within GUI
		logscaling  = boolean, if True (default) the logarithm of input_image is displayed
		colormap    = matplotlib color scheme used in the display
		interpolation = matplotlib interpolation scheme used for the display
		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# clear existing figure
		# plt.clf()

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# save input image for later use
		self.roi_obj.input_image = input_image

		ax = plt.subplot(111) 
		plt.subplots_adjust(left=0.05, bottom=0.2)

		# print out some instructions
		plt.suptitle('Use the slider bar to select ROIs, close the plotting window when satisfied.')

		# initial threshold value
		thres0 = 0.0

		# create a figure object
		figure_obj = plt.imshow(input_image,interpolation=interpolation)
		figure_obj.set_cmap(colormap)

		# prepare the slider bar
		thresxcolor = 'lightgoldenrodyellow'
		thresxamp  = plt.axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
		maxthreshold=np.floor(np.amax(input_image)) # maximum of slider
		sthres = plt.Slider(thresxamp, 'Threshold', 0.0, maxthreshold, valinit=thres0)

		textBox=plt.figtext(0.50, 0.065, 'Multiplier: 1.0',verticalalignment='center')

		# define what happens when the slider is touched
		def update(val):
			# parse a threshold from the slider
			thres     = sthres.val*thresMultiplier.factor
			# median filter the image
			newmatrix = signal.medfilt2d(input_image, kernel_size=kernel_size)
			# set pixels below the threshold to zero
			belowthres_indices = newmatrix < thres
			newmatrix[belowthres_indices] = 0
			# identify connected regions (this is already the roi_matrix)
			self.roi_obj.roi_matrix,numfoundrois = measurements.label(newmatrix)
			print str(numfoundrois) + ' ROIs found!'	
			figure_obj.set_data(newmatrix)
			plt.draw()

		# Buttons for changing multiplier for the value of slider
		class thresMultiplierClass:
			factor = 1.0;
			def __new__(cls):
				return self.factor
			def increase(self,event):
				self.factor *=2.0
				textBox.set_text('Multiplier: ' + str(self.factor))
				return self.factor
			def decrease(self,event):
				self.factor /=2.0
				textBox.set_text('Multiplier: ' + str(self.factor))
				return self.factor

		# call the update function when the slider is touched
		sthres.on_changed(update)

		thresMultiplier = thresMultiplierClass()
		axincrease   = plt.axes([0.8, 0.05, 0.05, 0.03])
		axdecrease   = plt.axes([0.7, 0.05, 0.05, 0.03])
		bnincrease    = Button(axincrease, 'x 2')
		bndecrease    = Button(axdecrease, '/ 2')
		bnincrease.on_clicked(thresMultiplier.increase)  # First change threshold
		bnincrease.on_clicked(update)		         # Then update image
		bndecrease.on_clicked(thresMultiplier.decrease)
		bndecrease.on_clicked(update)
		# ADDITION ENDS

		plt.show()

		# assign the defined rois to the roi_object class
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = xrs_rois.convert_matrix_rois_to_inds(self.roi_obj.roi_matrix)
		self.roi_obj.kind           = 'auto'
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(self.roi_obj.indices)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(self.roi_obj.indices)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def get_polygon_rois(self,input_image,modind=-1,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by clicking arbitrary number of points on a 2D image:
		LEFT CLICK to define the corner points of polygon, 
		MIDDLE CLICK to finish current ROI and move to the next ROI,
		RIGHT CLICK to cancel the previous point of polygon
		input_object   = 2D array, scan_object, or dictionary of scans to define the ROIs from
		modind	     = integer to identify module, if -1 (default), no module info will be in title (the case of one big image)
		logscaling     = boolean, to determine wether the image is shown on a log-scale (default = True)

		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# check that the input is a 2d matrix
		if not len(input_image.shape) == 2:
			print 'Please provide a 2D numpy array as input!'
			return

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# save input image for later use
		self.roi_obj.input_image = input_image

		# prepare a figure
		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		
		moduleNames='VD:','HR:','VU:','HL:','VB:','HB:',''  # for plot title

		# Initialize suptitle, which will be updated
		titlestring = ''
		titleInst=plt.suptitle(titlestring)

		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		# generate an image to be displayed
		figure_obj = plt.imshow(input_image,interpolation=interpolation)

		# set the colormap for the image
		figure_obj.set_cmap(colormap)

		rois   = []
		class Index:
			ind  = 1
			def next(self, event):

				titlestring = '%s next ROI is Nr. %02d:\n Left button to new points, middle to finish ROI. Hit \'Finish\' to end with this image.' % (moduleNames[modind], self.ind)
				titleInst.set_text(titlestring) # Update title

				# Try needed, as FINISH button closes the figure and ginput() generates _tkinter.TclError 
				try:
					one_roi = define_polygon_roi(input_image.shape)
					for index in one_roi:
						input_image[index[0],index[1]] += 1.0e6
					figure_obj.set_data(input_image)
					plt.hold(True)
					plt.draw()
					rois.append(one_roi)
					self.ind += 1
					# Begin defining the next ROI right after current
					self.next(self)
				except KeyboardInterrupt:	# to prevent "dead" figures
					plt.close()
					pass		
				except:
					pass
			def prev(self, event):
				self.ind -= 1
				for index in rois[-1]:
					input_image[index[0],index[1]] -= 1.0e6
				figure_obj.set_data(input_image)
				plt.hold(True)
				plt.draw()
				rois.pop()
				self.next(self)

			def close(self, event):
				plt.close()
			def dmy(self, event):
				pass # adding a dummy function for the dummy button

		callback = Index()	
		axprev   = plt.axes([0.5, 0.05, 0.1, 0.075])
		axnext   = plt.axes([0.61, 0.05, 0.1, 0.075])
		axclose  = plt.axes([0.72, 0.05, 0.1, 0.075])
		axdmy    = plt.axes([0.001, 0.001, 0.001, 0.001]) # for some reason the first botton disappears when clicked
		bdmy     = Button(axdmy,'')                        # which is why I am including a dummy button here
		bdmy.on_clicked(callback.dmy)                   # this way, the real buttons work
		bnext    = Button(axnext, 'Next')
		bnext.on_clicked(callback.next)
		bprev    = Button(axprev, 'Back')
		bprev.on_clicked(callback.prev)
		bclose   = Button(axclose, 'Finish')
		bclose.on_clicked(callback.close)

		# START: initiate NEXT button press
		callback.next(self)		

		# assign the defined rois to the roi_object class
		self.roi_obj.roi_matrix     = xrs_rois.convert_inds_to_matrix(rois,input_image.shape)
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = rois 
		self.roi_obj.kind           = 'polygon'
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(rois)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(rois)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def show_rois(self,interpolation='nearest'):
		"""
		Creates a figure with the defined ROIs as numbered boxes on it.
		"""
		roi_matrix = self.roi_obj.roi_matrix

		# check if there are ROIs defined
		if not np.any(roi_matrix):
			print 'Please select some rois first.'

		# make a figure
		else:
			plt.cla()
			plt.imshow(roi_matrix,interpolation=interpolation)
			plt.xlabel('x-axis [pixel]')
			plt.ylabel('y-axis [pixel]')

			# add a label with the number to the center of each ROI
			for ii in range(int(np.amax(roi_matrix))):
				# find center of the ROI and write the label
				inds    = np.where(roi_matrix[:,:] == ii+1)
				xcenter = np.mean(inds[1])
				ycenter = np.mean(inds[0])
				string  = '%02d' % (ii+1)
				plt.text(xcenter,ycenter,string)

	def import_simo_style_rois(self,roiList,detImageShape=(512,768)):
		""" **import_simo_style_rois**
		Converts Simo-style ROIs to the conventions used here.

		Arguments:
		----------
		roiList (list): List of tuples that have [(xmin, xmax, ymin, ymax), (xmin, xmax, ymin, ymax), ...].
		detImageShape (tuple): Shape of the detector image (for convertion to roiMatrix)
		"""
		indices = []
		for roi in roiList:
			inds = []
			for ii in range(roi[0],roi[1]):
				for jj in range(roi[2],roi[3]):
					inds.append((jj,ii))
			indices.append(inds)
		# assign the defined rois to the roi_object class
		if detImageShape:
			self.roi_obj.roi_matrix     = xrs_rois.convert_inds_to_matrix(indices,detImageShape)
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = indices
		self.roi_obj.kind           = 'simoStyle'
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(indices)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(indices)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def refine_pw_rois(self,roi_obj,pw_data,n_components=2,method='nnma',cov_thresh=-1):
		"""
		**refine_pw_rois**

		Use decomposition of pixelwise data for each ROI to find which of the pixels holds
		data from the sample and which one only has background. 

		Args:
		-----
		roi_obj (xrs_rois.roi_object): ROI object to be refined
		pw_data (list): list containing one 2D numpy array per ROI holding pixel-wise signals
		n_components (int): number of components in the decomposition
		method (string): keyword describing which decomposition to be used ('pca', 'ica', 'nnma')
		"""
		# check if available method is used
		avail_methods = ['pca','ica','nnma']
		if not method in avail_methods:
			print('Please use one of the following methods: ' + str(avail_methods) + '!')
			return

		# check if scikit learn is available
		try:
			from sklearn.decomposition import FastICA, PCA, ProjectedGradientNMF
		except ImportError:
			raise ImportError('Please install the scikit-learn package to use this feature.')
			return

		counter = 0
		new_indices = []
		for data in pw_data: # go through each ROI

			# decompose data, choose method
			if method == 'nnma': # non negative matrix factorisation
				nnm = ProjectedGradientNMF(n_components=n_components)
				N   = nnm.fit_transform(data)
			elif method == 'pca': # principal component analysis
				pca = PCA(n_components=n_components)
				N = pca.fit_transform(data)
			elif method == 'ica': # independent component analysis
				ica = FastICA(n_components=n_components)
				N = ica.fit_transform(data)

			# let user decide which component belongs to the data:
			user_choise = 0
			plt.cla()
			plt.title('Click on component you are interested in.')
			legendstr = []
			for ii in range(n_components):
				plt.plot(N[:,ii])
				legendstr.append('Component No. %01d' %ii)
			plt.legend(legendstr)
			plt.xlabel('points along scan')
			plt.ylabel('intensity [arb. units]')
			user_input = np.array(plt.ginput(1,timeout=-1)[0])

			# which curve was chosen
			nearest_points = [(np.abs(N[:,ii]-user_input[1])).argmin() for ii in range(n_components)]
			user_choice = (np.abs(nearest_points-user_input[0])).argmin()

			# find covariance for all pixels with user choice
			covariance = np.array([])
			for ii in range(len(data[0,:])):
				covariance = np.append(covariance, np.cov(data[:,ii],N[:,user_choice])[0,0])

			# plot covariance, let user choose the the cutoff in y direction
			plt.cla()
			plt.title('Click to define a y-threshold.')
			plt.plot(covariance)
			plt.xlabel('pixels in ROI')
			plt.ylabel('covariance [arb. units]')
			if cov_thresh < 0:
				user_cutoff = np.array(plt.ginput(1,timeout=-1)[0])
			elif cov_thresh>0 and isinstance(cov_thresh,int):
				if len(covariance) < cov_thresh:
					print('ROI has fewer pixels than cov_thresh, will break here.')
					return
				else:
					user_cutoff = np.array([0.0, np.sort(covariance)[-cov_thresh]])
			else:
				print('Please provide cov_thresh as positive integer!')

			# find the ROI indices above the cutoff, reassign ROI indices
			inds = covariance >= user_cutoff[1]
			refined_indices = []
			for ii in range(len(roi_obj.indices[counter])):
				if inds[ii]:
					refined_indices.append(roi_obj.indices[counter][ii])
			new_indices.append(refined_indices)

			# end loop
			counter += 1

		# reassign ROI object
		self.roi_obj.roi_matrix     = xrs_rois.convert_inds_to_matrix(new_indices,self.roi_obj.input_image.shape)
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = new_indices 
		self.roi_obj.kind           = 'refined'
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(new_indices)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(new_indices)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)


def define_lin_roi(height,image_shape,verbose=False):
	"""
	Lets you pick 2 points on a current image and returns a linear ROI of
	height (2*height+1).
	height      = number of pixels that define the height of the ROI
	image_shape = tuple with shape of the current image (i.e. (256,256))
	"""
	endpoints = list(np.round(plt.ginput(2,timeout=-1)))

	# check that selected points are in the image
	for point in endpoints:
		if point[0] < 0.0:
			point = 0
		if point[0] > image_shape[1]:
			point = image_shape[1]
		if point[1] < 0.0:
			point[1] = 0
		if point[1] > image_shape[0]:
			point[1] = image_shape[0]

	# check that point 2 is bigger than point 1 in the x direction
	if endpoints[1][0]< endpoints[0][0]:
		endpoints[0],endpoints[1] = endpoints[1],endpoints[0]

	# print the limits of the rectangle in the shell
	if not verbose:
		print 'The selected points are: ', [[endpoints[0][1],endpoints[0][0]],[endpoints[1][1],endpoints[1][0]]]

	roix = np.arange(endpoints[0][0],endpoints[1][0])
	roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]
	roi  = []
	height = np.arange(-height,height)
	for n in range(len(roix)):
		for ind in height:
			roi.append((roiy[n]+ind,roix[n]))
	return roi

def define_zoom_roi(input_image,verbose=False):
	"""
	Parses the current figure limits and uses them to define a rectangle of "
	roi_number"s in the matrix given by roi_matrix.
	input_image  = unzoomed figure
	roi_matrix   = current roi matrix which will be altered 
	"""
	# input_image.shape prints  (y-length, x-length)
	# parse the figure limits from the current zoom
	# limits = np.round(plt.axis()) # x-min, x-max, y-max, y-min
	frac_limits = plt.axis() # x-min, x-max, y-max, y-min as floats
	limits = np.array([np.ceil(frac_limits[0]), np.floor(frac_limits[1]+0.5), np.floor(frac_limits[2]+0.5), np.ceil(frac_limits[3])]) #  x-min, x-max, y-max, y-min as ints
	# check that selected zoom area is not outside the image
	inds = limits < 0
	limits[inds] = 0
	if limits[1] > input_image.shape[1]:
		limits[1] = input_image.shape[1]
	if limits[2] > input_image.shape[0]:
		limits[2] = input_image.shape[0]

	# sort the limits in ascenging order
	limitsy = limits[2:4] # vertical
	limitsy.sort()
	limitsx = limits[0:2] # horizontal
	limitsx.sort()

	# print the limits of the rectangle in the shell
	if not verbose:
		print 'The selected limits are: ', limitsx, limitsy

	# prepare a n*m matrix with one roi
	T = np.zeros(input_image.shape)
	T[limitsy[0]:limitsy[1],limitsx[0]:limitsx[1]] = 1
	indsy,indsx = np.where(T == 1)

	roi  = []
	for n in range(len(indsx)):
		roi.append((indsy[n],indsx[n]))
	return roi

def show_rois(roi_matrix):
        """
        Creates a figure with the defined ROIs as numbered boxes on it.
        """
        roi_matrix 

        # check if there are ROIs defined
        if not np.any(roi_matrix):
                print 'Please select some rois first.'

        # make a figure
        else:
                plt.ioff()
                plt.imshow(roi_matrix)
                plt.xlabel('x-axis [pixel]')
                plt.ylabel('y-axis [pixel]')

                # add a label with the number to the center of each ROI
                for ii in range(int(np.amax(roi_matrix))):
                        # find center of the ROI and write the label
                        inds    = np.where(roi_matrix[:,:] == ii+1)
                        xcenter = np.mean(inds[1])
                        ycenter = np.mean(inds[0])
                        string  = '%02d' % (ii+1)
                        plt.text(xcenter,ycenter,string)
                plt.show()
                plt.ion()



def test_roifinder(roi_type_str, imagesize = [512,768], scan = None ):
    """
    Runs the roi_finder class on a random image of given type for 
    testing purposes.
    scan[0] = absolute path to a spec file
    scan[1] = energycolumn='energy'
    scan[2] = monitorcolumn='kap4dio'
    scan[3] = scan number from which to take images
    """
    strings = ['zoom','linear','auto']
    if not roi_type_str in strings:
	    print 'Only ' + str(strings) + ' testable, choose one of them!' 
	    return

    # make a random image 
    if not scan:
	    rand_image = np.random.rand(imagesize[0],imagesize[1])
    else: 
	    import xrs_read, xrs_utilities
	    read_obj = xrs_read.read_id20(scan[0],energycolumn=scan[1],monitorcolumn=scan[2])
	    read_obj.loadelastic(scan[3])
	    key = 'Scan%03d' % scan[3]
	    rand_image = xrs_utilities.sumx(read_obj.scans[key].edfmats)

    # create a roi_finder object 
    roi_finder_obj = roi_finder()

    if roi_type_str == 'zoom':
	    roi_finder_obj.get_zoom_rois(rand_image,logscaling=True,colormap='jet',interpolation='nearest')
	    roi_finder_obj.show_rois()

    elif roi_type_str == 'linear':
	    roi_finder_obj.get_linear_rois(rand_image,logscaling=True,height=5,colormap='jet',interpolation='nearest')
	    roi_finder_obj.show_rois()
	    
    elif roi_type_str == 'auto':
	    roi_finder_obj.get_auto_rois(rand_image,kernel_size=5,threshold=1.0,logscaling=True,colormap='jet',interpolation='nearest')
	    roi_finder_obj.show_rois()



def create_diff_image(scans,scannumbers,energy_keV):
	"""
	Returns a summed image from all scans with numbers 'scannumbers'.
	scans       = dictionary of objects from the scan-class
	scannumbers = single scannumber, or list of scannumbers from which an image should be constructed
	"""
	# make 'scannumbers' iterable (even if it is just an integer)
	numbers = []
	if not isinstance(scannumbers,list):
		numbers.append(scannumbers)
	else:
		numbers = scannumbers

	key = 'Scan%03d' % numbers[0]
	below_image = np.zeros_like(scans[key].edfmats[0,:,:])
	above_image = np.zeros_like(scans[key].edfmats[0,:,:])

	# find indices below and above 'energy'
	below_inds = scans[key].energy < energy_keV
	above_inds = scans[key].energy > energy_keV
	for number in numbers:
		key = 'Scan%03d' % number
		for ii in below_inds:
			below_image += scans[key].edfmats[ii,:,:]
		for ii in above_inds:
			above_image += scans[key].edfmats[ii,:,:]

	return (above_image - below_image)





def define_polygon_roi(image_shape,verbose=False):
	"""
	Define a polygon shaped ROI from a current image by 
	selecting points.
	"""

	tmptuple=plt.ginput(0,timeout=-1,show_clicks=True)
	pnts=np.array(tmptuple)
	if verbose:
		print 'The selected points are:' 
		print pnts
	# Create the polygon path from points clicked
	path=Path(pnts)

	# bounding box: everything outside this is necessarily zero
	bbx = np.arange( np.floor(min( pnts[:,0])), np.ceil(max( pnts[:,0])))
	bby = np.arange( np.floor(min( pnts[:,1])), np.ceil(max( pnts[:,1])))
	
	roi=[]
	for ii in bbx:
		for jj in bby:
			# Test point
			if path.contains_point([ii, jj]):
				roi.append( (jj, ii))  

	return roi

