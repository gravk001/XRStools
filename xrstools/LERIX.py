from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

try:
    import os, sys, re, time, h5py, glob
    import numpy as np
    from numpy import array
    import pandas as pd
    from  dateutil.parser import parse as dateparse
    from datetime import datetime
except:
    print('LERIX requires the following modules: os,sys,time,h5py,glob,numpy,pandas.dateutil.parser,datetime')

try:
    from XRStools import xrs_utilities, xrs_scans #xrs_scans used to make a scan group for XRStools and utilities used to find centre of gaussians
except:
    print("XRStools is required to run LERIX! Please see: 'https://github.com/christophsahle/XRStools.git'")

TINY = 1.e-7
MAX_FILESIZE = 100*1024*1024  # 100 Mb limit
COMMENTCHARS = '#;%*!$'
NAME_MATCH = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)*$").match
VALID_SNAME_CHARS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
VALID_NAME_CHARS = '.%s' % VALID_SNAME_CHARS
RESERVED_WORDS = ('and', 'as', 'assert', 'break', 'class', 'continue',
                'def', 'del', 'elif', 'else', 'eval', 'except', 'exec',
                'execfile', 'finally', 'for', 'from', 'global', 'if',
                'import', 'in', 'is', 'lambda', 'not', 'or', 'pass',
                'print', 'raise', 'return', 'try', 'while', 'with',
                'group', 'end', 'endwhile', 'endif', 'endfor', 'endtry',
                'enddef', 'True', 'False', 'None')

################################################################################
# Functions to get the header attributes
################################################################################
class Lerix:

    def __init__(self):
        self.scans         = {}
        self.data          = {} #np dictionary of arrays separated into their column headers
        self.header_attrs  = {} #dictionary of useful information from scan files, inc. e0, comments, scan_time+date
        self.key           = {'Analyzer1':0, 'Analyzer2':1, 'Analyzer3':2,'Analyzer4':3,'Analyzer5':4,
                            'Analyzer6':5,'Analyzer7':6,'Analyzer8':7,'Analyzer9':8,'Analyzer10':9,'Analyzer11':10,'Analyzer12':11
                            ,'Analyzer13':12,'Analyzer14':13,'Analyzer15':14,'Analyzer16':15,'Analyzer17':16,'Analyzer18':17,'Analyzer19':18}
        self.elastic_scans = []
        self.elastic_name  = 'elastic'
        self.nixs_scans    = []
        self.nixs_name     = 'nixs'
        self.wide_scans    = []
        self.wide_name     = 'wide'
        self.scan_name     = []
        self.sample_name   = []
        self.eloss_avg     = [] #elastic eloss average
        self.signals_avg   = [] #elastic signals average used to plot analyzer resolutions at the end
        self.energy        = []
        self.signals       = []
        self.errors        = []
        self.is_checked    = [] #inserted later to save a list of the chosen analyzers after using .plot_data() save function
        self.tth           = []
        self.resolution    = {}
        self.E0            = []
        self.cenom         = []
        self.cenom_dict    = {}

    ################################################################################
    # Get Ascii Info - parse a file and return the key details
    ################################################################################
    def getfloats(self, txt, allow_times=True):
        """
        function goes through a line and returns the line as a list of strings
        """
        words = [w.strip() for w in txt.replace(',', ' ').split()]
        mktime = time.mktime
        for i, w in enumerate(words):
            val = None
            try:
                val = float(w)
            except ValueError:
                try:
                    val = mktime(dateparse(w).timetuple())
                except ValueError:
                    pass
        words[i] = val
        return(words)

    def colname(self, txt):
        """Function to replace bad characters with '_''s making a line of strings
        easier to handle."""
        return self.fixName(txt.strip().lower()).replace('.', '_')

    def isValidName(self, filename):
        """Function checks that a filename isn't in the list of reserved pythonic
        words. Returns corrected name or False"""
        if filename in RESERVED_WORDS:
            return False
        tnam = filename[:].lower()
        return NAME_MATCH(tnam) is not None

    def fixName(self, filename, allow_dot=True):
        if self.isValidName(filename):
            return filename
        if self.isValidName('_%s' % filename):
            return '_%s' % filename
        chars = []
        valid_chars = VALID_SNAME_CHARS
        if allow_dot:
            valid_chars = VALID_NAME_CHARS
        for s in filename:
            if s not in valid_chars:
                s = '_'
            chars.append(s)
        filename = ''.join(chars)
        # last check (name may begin with a number or .)
        if not self.isValidName(filename):
            filename = '_%s' % filename
        return filename

    def strip_headers(self, headers):
        #reorganise the headers and remove superfluous lines and commentchars
        header = []
        for hline in headers:
            hline = hline.strip().replace('\t', ' ')
            if len(hline) < 1:
                continue
            if hline[0] in COMMENTCHARS:
                hline = hline[1:].lstrip() #assumes reading l2r
            if len(hline) <1:
                continue
            header.append(hline)
        return(header)

    def separate_infile(self, text):
        """Function parses an 20ID ASCII file in reverse and separates it into
        headers, footers and data"""
        _labelline = None
        ncol = None
        dat, footers, headers = [], [], []
        try:
            text.reverse()
        except:
            text[::-1]
        section = 'FOOTER'
        for line in text:
            line = line.strip()
            if len(line) < 1: #remove any blank lines
                continue
            if section == 'FOOTER' and not None in self.getfloats(line):
                section = 'DATA'
            elif section == 'DATA' and None in self.getfloats(line):
                section = 'HEADER'
                _labelline = line
                if _labelline[0] in COMMENTCHARS:
                    _labelline = _labelline[1:].strip()
            if section == 'FOOTER': #reading footers but not using them currently
                footers.append(line)
            elif section == 'HEADER':
                headers.append(line)
            elif section == 'DATA':
                rowdat  = self.getfloats(line)
                if ncol is None:
                    ncol = len(rowdat)
                if ncol == len(rowdat):
                    dat.append(rowdat)
        return(headers, dat, footers)

    def pull_id20attrs(self, header):
        """Function takes headers of 20ID ASCII and parses it for key information
        for header_attrs - N.B. could be shortened by looping through a list of
        important key words rather than doing one by one."""
        bounds, steps, int_times = [], [], []
        header_attrs = {}
        line = -2
        #iterate through the header and pull out useful information and send it to header_attrs Dictionary
        for hhline in map(str.lower,header):
            line = line + 1 #counting to return the user comments which are on the next line
            try:
                if str(header[comment_line].strip()) == 'Scan config:':
                    header_attrs['User Comments'] = ""
                    pass
                else:
                    header_attrs['User Comments'] = str(header[comment_line].strip())
            except:
                pass
            if hhline.startswith('beamline'):
                words = hhline.split('beamline',1)
                header_attrs['beamline'] = str(words[1].strip())
            elif hhline.startswith('e0'):
                if ':' in hhline:
                    words = hhline.split(':',1)
                    header_attrs[words[0]] = float(words[1].strip(' ').split(' ',1)[0])
                elif '=' in hhline:
                    words = hhline.split('=',1)
                    header_attrs[words[0]] = float(words[1].strip(' ').split(' ',1)[0])
            elif hhline.startswith('user comment'):
                comment_line = line
            elif "scan time" in hhline:
                #search for scan date and time see: https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
                try:
                    words = hhline.split('scan time',1)
                    header_attrs['scan_time'] = datetime.strptime(words[1].strip(), '%H hrs %M min %S sec.').time()
                    header_attrs['scan_date'] = datetime.strptime(words[0].split('panel',1)[1].strip().strip(';'), '%m/%d/%Y  %I:%M:%S %p').date()
                except:
                    continue
            elif "scan bounds" in hhline:
                words = hhline.split('scan bounds',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        bounds.append(float(i))
                    except:
                        pass
                header_attrs['scan_bounds'] = bounds
            elif "scan step(s)" in hhline:
                words = hhline.split('scan step(s)',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        steps.append(float(i))
                    except:
                        pass
                header_attrs['scan_steps'] = steps
            elif "integration times" in hhline:
                words = hhline.split('integration times',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        int_times.append(float(i))
                    except:
                        pass
                header_attrs['int_times'] = int_times
        return(header_attrs)

    def get_col_headers(self, header):
        col_headers = []
        for i in self.colname(header[0]).split('___'): #need three _ to work
            if not i:
                continue
            col_headers.append(i.strip('_'))
        return(col_headers)

    def scan_info(self, f):
        """get the scan number, name, type and file extention from the title of
        the scan assuming typical format e.g. elastic.0001, nixs.0001"""
        f = os.path.basename(f) #this allows both directories and files to be passed to get scan_info
        fn,fext = os.path.splitext(f)
        if str.lower(fn)==str.lower(self.nixs_name):
            scan_type = 'nixs'
        elif str.lower(fn)==str.lower(self.elastic_name):
            scan_type = 'elastic'
        elif str.lower(fn)==str.lower(self.wide_name):
            scan_type = 'wide'
        else:
            print(""">> LERIX >> WARNING  \n You have probably called the scan_info
            function without specifying a correct \n <class>.nixs/wide/elastic_name if you
            are calling scan_info manually - you can change this by setting:\n\n
            <class>.nixs_name = '<nixs_name>'""")
            sys.exit()
        scan_number = fext.lstrip('.')
        scan_number = int(scan_number)
        scan_name = scan_type + '%04d' %scan_number
        return(scan_number, scan_name, scan_type, f)

    def sort_dir(self, dir):
        """Returns a list of directory contents after filtering out scans without
        the correct format or size e.g. 'elastic.0001, nixs.0001 '"""
        dir_scans = []
        for file in os.listdir(dir):
            file_lc = str.lower(file)
            fn,fext = os.path.splitext(file_lc)
            if not file_lc.startswith('.'):
                    if fext.lstrip('.').isdigit():
                        if not os.stat(dir + '/' + file).st_size > 8000:
                            print("{} {}".format(">> >> Warning!! skipped empty scan (<8KB): ", file))
                            continue
                        elif not os.stat(dir + '/' + file).st_size < MAX_FILESIZE:
                            print("{} {}".format(">> >> Warning!! skipped huge scan (>100MB): ", file))
                            continue
                        else:
                            if fn==self.nixs_name:
                                dir_scans.append(file)
                            elif fn==self.elastic_name:
                                dir_scans.append(file)
                            elif fn==self.wide_name:
                                dir_scans.append(file)
        sorted_dir = sorted(dir_scans, key=lambda x: os.path.splitext(x)[1])
        return sorted_dir

    def isValidDir(self,dir):
        """Show that the scan directory is valid, that the directory holds a scan
        with the correct elastic name, nixs name and then let the user know if it
        has not found a wide scan. Returns True if valid directory."""
        if not os.path.isdir(dir):
            print('Check the directory you have supplied')
            return False
        elif not os.path.isfile(dir+'/'+self.elastic_name+'.0001'):
            print("The directory you supplied does not have a elastic.0001 file!!! \n If your elastic scan has a different name, please specify as: 'elastic_name'")
            return False
        elif not os.path.isfile(dir+'/'+self.nixs_name+'.0001'):
            print("The directory you supplied does not have a NIXS.0001 file!!! \n If your raman scan has a different name, please specify as: 'NIXS_name'")
            return False
        elif not os.path.isfile(dir+'/'+self.wide_name+'.0001'):
            print("No wide scans found. Continuing...")
            return True
        else:
            return True

    def plot_data(self,analyzer=False):
        """<classObj>.plot_data() Function that can be called to plot the eloss
        data for each channel and build an average by clicking a button.
        Does not require matplotlib >2.1"""
        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, Button, Cursor
        channels = []
        for analyzer in self.resolution:
            if analyzer.startswith('Analyzer'):
                if self.resolution[analyzer] < 1.0:
                    channels.append(int(analyzer.lstrip('Analyzer'))-1)
        data = np.average(self.signals[:,channels],axis=1)
        fig, ax = plt.subplots()
        ax.plot(self.eloss, data, lw=2)
        ax.set_xlabel('Energy Loss (eV)')
        ax.set_ylabel('S(q,w) [1/eV]')
        ax.set_title('Plotting Raman Analysers')
        plt.subplots_adjust(left=0.3)
        checkbuttonaxis = plt.axes([0.02, 0.15, 0.2, 0.8])
        anlabels, anvals = list(self.key), (False,)*len(list(self.key))
        anstates = dict(zip(anlabels,anvals))
        analyzers = CheckButtons(checkbuttonaxis, anlabels, anvals)
        buttonaxis = plt.axes([0.01, 0.01, 0.3, 0.09])
        bStatus  = Button(buttonaxis,'Save Averaged Analyzers')

        def onclick(label):
            """
            Tell the user what they have clicked - also good for de-bugging
            """
            anstates[label] = not anstates[label]
            print('un'*(not anstates[label]) + 'checked %s' %label)
            func()

        def savebutton(val):
            from PyQt4.QtGui import QApplication, QWidget, QFileDialog
            if not self.is_checked:
                print('please select your chosen analysers first!')
            else:
                print('selected analysers (python counting):  ', self.is_checked)
                save_signals = np.average(self.signals[:,self.is_checked],axis=1)
                save_errors = np.average(self.errors[:,self.is_checked],axis=1)
                df = pd.DataFrame(list(zip(self.eloss,save_signals,save_errors)), columns=['eloss','signals','errors'])
                print(df)
                try:
                    w = QWidget() #bug 'start a Qapp before a QPaintDevice py36 mac'
                    filename = str(QFileDialog.getSaveFileName(w, 'Save Analyzer Average','Result.csv'))
                    df.to_csv(filename,sep=',',na_rep='nan')
                    print('Saved as: ',filename)
                except:
                    print("{} {}".format(">> Warning >>", "file save was unsuccessful"))

        def func():
            ax.clear()
            self.is_checked = []
            for ii in anlabels:
                if anstates[ii]:
                    self.is_checked.append(self.key[ii])
            ax.plot(self.eloss, np.average(self.signals[:,self.is_checked],axis=1),lw=2)
            cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
            ax.autoscale(True)
            ax.set_xlabel('Energy Loss (eV)')
            ax.set_title('Plotting Raman Analysers')
            plt.draw()

        bStatus.on_clicked(savebutton)
        analyzers.on_clicked(onclick)
        cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
        plt.show()

    def write_H5scanData(self,dir,H5file,H5name,averaged='False'):
        """Writes all the scan information into a H5 file named after the sample name. inside
        this H5 directory scans are split into elastic and NIXS and the averaged scans. No support
        yet for wide scans"""
        g = H5file.create_group(H5name) #H5 subgroup with the name of the sample
        H5_ela = g.create_group('elastic') #H5 subgroup for elastics
        H5_xrs = g.create_group('XRS')     #H5 subgroup for NIXS
        all_scans = self.elastic_scans+self.nixs_scans
        for file in all_scans:
            scan_info = self.scan_info(file)
            if scan_info[2] == 'elastic':
                h5group = H5_ela.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                h5group.create_dataset("cenoms",data=self.scans[scan_info[1]].cenom)
            elif scan_info[2]=='nixs':
                h5group = H5_xrs.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("eloss",data=self.scans[scan_info[1]].eloss)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                h5group.create_dataset("tth",data=self.scans[scan_info[1]].tth)

        g.create_dataset("energy",data=self.energy)
        g.create_dataset("signals",data=self.signals)
        g.create_dataset("eloss",data=self.eloss)
        g.create_dataset("errors",data=self.errors)
        g.create_dataset("tth",data=self.tth)
        g.create_dataset("Mean Resolutions", data=np.array(self.resolution.items()))

        #Never forget to close an open H5 file!!!
        H5file.close()

    ################################################################################
    # Read Scan
    ################################################################################
    def get_cenoms(self, file):
        """Internal Function to get the centre of mass of the elastic peak and
        the E0 for each elastic scan using XRStools"""
        cenom_list, cenom_analyzers = [], []
        scan_info = self.scan_info(file)
        for analyzer in range(19): #The analyzer channels in the scan ASCII
            #self.scans[scan_info[1]].cenom.append(xrs_utilities.find_center_of_mass(self.scans[scan_info[1]].energy,self.scans[scan_info[1]].signals[:,analyzer]))
            cenom_analyzers.append(xrs_utilities.find_center_of_mass(self.scans[scan_info[1]].energy,self.scans[scan_info[1]].signals[:,analyzer]))
            self.scans[scan_info[1]].cenom = np.mean(cenom_analyzers) #list of cenoms not necessary here since energy is the same for each analyser in a scan.
        cenom_list.append(self.scans[scan_info[1]].cenom)
        # self.cenom = [sum(a)/len(a) for a in zip(*cenom_list)]
        print(cenom_list)
        self.cenom = np.mean(cenom_list)
        self.E0 = np.mean(self.cenom)/1e3

    def readscan_20ID(self, file, valid_elastic=False):
        """Read an ID20-type ASCII file and return header attributes and data as
        a dictionary. Takes a file path.
        header_attrs -> int_times, scan_steps, scan_bounds, e0, comments, beamline,
                        scan_time, scan_date
        data         -> dictionary of np.array (float64) with callable column names"""
        scan_info = self.scan_info(file)
        qixs_list = []
        f = open(file, "r") #read starts here
        text = f.read()
        text = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
        headers, dat, footers = self.separate_infile(text)
        try:
            dat = [map(list,zip(*dat))[i][::-1] for i in range(len(dat[1]))] # this function does the inverse ([::-1]) transposition of the dat object, doesn't seem to work in windows
        except:
            dat = [list(map(list,zip(*dat)))[i][::-1] for i in range(len(dat[1]))]
        names = self.get_col_headers(self.strip_headers(headers)) #returns a list of names in the order found in the data file.
        data = pd.DataFrame(np.array(dat).T,columns = np.array(names).T, dtype='float64') #returns a pandas array with the data arranged into labelled columns
        for column in sorted(data.columns): #sort by name so that analyzers are in correct (numerical) order
            if not column.rfind('i0') == -1:
                self.scans[scan_info[1]].monitor = data[column].values
            if not column.rfind('__alt') == -1:
                self.scans[scan_info[1]].energy = data[column].values
            if not column.rfind('qixs') == -1:
                qixs_list.append(column)
        self.scans[scan_info[1]].signals = data[qixs_list].values
        self.scans[scan_info[1]].errors  = np.sqrt(np.absolute(self.scans[scan_info[1]].signals))
        if scan_info[2]=='elastic':
            self.get_cenoms(file)
            self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,self.scans[scan_info[1]].cenom)
        elif scan_info[2]=='nixs' or scan_info[2]=='wide':
            #create empty array with shape energy.v.signals
            eloss = np.zeros(self.scans[scan_info[1]].signals.shape)
            self.scans[scan_info[1]].tth = list(range(9,180,9)) #assign tth to each scan
            self.tth = list(range(9,180,9)) #assign tth to self
            if valid_elastic:
                print('>>>>>>> VALID ELASTIC')
                self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,self.scans['elastic%04d'%scan_info[0]].cenom)
            elif not valid_elastic:
                print('>>>>>>> NO VALID ELASTIC')
                try:
                    self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,self.cenom)
                except:
                    print('>> LERIX >> Reading a NIXS scan without any elastic scans reduces confidence in the result. LERIX is now taking E0 from the scan header file')
                    scan_attrs = self.pull_id20attrs(self.strip_headers(headers)) #get scan_attrs
                    self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,scan_attrs['e0'])
        f.close()

    def average_scans(self,scan_numbers='all'):
        """Function to calculate the average eloss, energy, signals and errors over
        all the read scans (default) or over a list of scans e.g. [1,3,5]"""
        energy_running = []
        signals_running = []
        eloss_running = []
        errors_running = []
        if scan_numbers=='all':
            for file in self.nixs_scans:
                scan_info = self.scan_info(file)
                energy_running.append(self.scans[scan_info[1]].energy)
                signals_running.append(self.scans[scan_info[1]].signals)
                eloss_running.append(self.scans[scan_info[1]].eloss)
                errors_running.append(self.scans[scan_info[1]].errors)
            self.energy = np.array([sum(a)/len(a) for a in zip(*energy_running)])
            self.signals = np.array([sum(a)/len(a) for a in zip(*signals_running)])
            self.eloss = np.array([sum(a)/len(a) for a in zip(*eloss_running)])
            self.errors = np.array([sum(a)/len(a) for a in zip(*errors_running)])
        elif type(scan_numbers) is list:
            scan_numbers[:] = [x - 1 for x in scan_numbers] #scan 1 will be the 0th item in the list
            chosen_scans = []
            for number in scan_numbers:
                scan_info = self.scan_info(self.nixs_scans[number])
                chosen_scans.append(scan_info[1])
            print("{} {}".format("Averaging scan numbers: ", chosen_scans))
            for scan in chosen_scans:
                energy_running.append(self.scans[scan_info[1]].energy)
                signals_running.append(self.scans[scan_info[1]].signals)
                eloss_running.append(self.scans[scan_info[1]].eloss)
                errors_running.append(self.scans[scan_info[1]].errors)
            self.energy = np.array([sum(a)/len(a) for a in zip(*energy_running)])
            self.signals = np.array([sum(a)/len(a) for a in zip(*signals_running)])
            self.eloss = np.array([sum(a)/len(a) for a in zip(*eloss_running)])
            self.errors = np.array([sum(a)/len(a) for a in zip(*errors_running)])
        else:
            print("scan_numbers must be blank, 'all' or a list of scan numbers e.g.[1,3,5]")
            sys.exit()

    def get_resolutions(self,scan_numbers):
        """Internal function to get the average resolution of each analyzer and
        average to give a self.resolution over the 19 analyzers. Returns a Dictionary
        of resolutions the mean and each analyser"""
        eloss_running_elastic = []
        signals_running_elastic = []
        if scan_numbers=='all':
            chosen_scans = []
            for number in range(len(self.elastic_scans)):
                scan_info = self.scan_info(self.elastic_scans[number])
                chosen_scans.append(scan_info[3])
        elif type(scan_numbers) is list:
            scan_numbers[:] = [x - 1 for x in scan_numbers] #scan 1 will be the 0th item in the list
            chosen_scans = []
            for number in scan_numbers:
                scan_info = self.scan_info(self.elastic_scans[number])
                chosen_scans.append(scan_info[3])
        else:
            print('scan numbers must be a list of the scans with correct length')
            return
        #populate lists with eloss and signals and then find the average over the whole range
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            eloss_running_elastic.append(self.scans[scan_info[1]].eloss)
            signals_running_elastic.append(self.scans[scan_info[1]].signals)
        self.eloss_avg = np.array([sum(a)/len(a) for a in zip(*eloss_running_elastic)])
        self.signals_avg = np.array([sum(a)/len(a) for a in zip(*signals_running_elastic)])
        #take these average values and find the average FWHM for each analyzer and then find the total average
        for file in chosen_scans:
            resolution = []
            skipped = []
            for analyzer in range(19
                try:
                    resolution.append(xrs_utilities.fwhm(self.eloss_avg, self.signals_avg[:,analyzer])[0])
                    self.resolution['Analyzer%s'%analyzer] = resolution[analyzer]
                except:
                    skipped.append(analyzer+1)
                    continue
        if len(skipped) > 1:
            print("{} {}".format("Skipped resolution for analyzer/s: ", list(set(skipped))))
        self.resolution['Resolution'] = round(np.mean(resolution),3)

    ################################################################################
    # Begin the reading
    ################################################################################
    def load_experiment(self,dir,nixs_name='NIXS',wide_name='wide',elastic_name='elastic',scan_numbers='all',H5=False,H5path=None,sample_name=None):
        """Function to load scan data from a typical APS 20ID Non-Resonant inelastic
        X-ray scattering experiment. With data in the form of elastic.0001, allign.0001
        and NIXS.0001. Function reteurns the averaged energy loss, signals, errors, E0
        and 2theta angles for the scans in the chosen directory."""

        #make sure the user inputs are all lower case for easy reading
        self.nixs_name = str.lower(nixs_name)
        self.wide_name = str.lower(wide_name)
        self.elastic_name = str.lower(elastic_name)

        #check dir location
        if not self.isValidDir(dir):
            print("{} {}".format(">> >> WARNING: ", "IO Error - check the directory name you have given"))
            sys.exit()
        else:
            pass

        #sort the directory so that scans are in order, determine number of scans
        #open list to be filled with the elastic/nixs scan names
        # number_of_scans = len(glob.glob(dir+'/'+self.nixs_name+'*'))-1 #removed because not called?!
        self.elastic_scans = []
        self.nixs_scans = []
        self.wide_scans = []
        #self.keys = {"eloss":np.array, "energy":np.array, "signals":np.array, "errors":np.array,"E0":np.float, "tth":np.array} #,"resolution":array }

        #split scans into NIXS and elastic and begin instance of XRStools scan class for each scan
        for file in self.sort_dir(dir):
                scan_info = self.scan_info(file)
                scan = xrs_scans.Scan()

                if scan_info[2]=='elastic':
                    self.elastic_scans.append(file)
                    self.scans[scan_info[1]] = scan #self.scans {} in class _init_
                    self.scans[scan_info[1]].scan_type = scan_info[2]
                    self.scans[scan_info[1]].scan_number = scan_info[0]

                if scan_info[2]=='nixs':
                    self.nixs_scans.append(file)
                    self.scans[scan_info[1]] = scan #self.scans {} in class _init_
                    self.scans[scan_info[1]].scan_type = scan_info[2]
                    self.scans[scan_info[1]].scan_number = scan_info[0]

                if scan_info[2]=='wide':
                    self.wide_scans.append(file)
                    self.scans[scan_info[1]] = scan #self.scans {} in class _init_
                    self.scans[scan_info[1]].scan_type = scan_info[2]
                    self.scans[scan_info[1]].scan_number = scan_info[0]

                else:
                    continue

        #read elastic scans first to calculate cenom
        # reading all the elastics in the file to improve E0 accuracy and get a
        # good grasp on the scan resolution
        for file in self.elastic_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading elastic scan: ", file))
            self.readscan_20ID(dir + '/' + file)

        print('I always read all the Elastic scans to improve Resolution and E0 Accuracy\n >> Type <class>.resolution to see the analyzer resolutions.')

        #Read NIXS scans - if there isn't a corresponing elastic scan, subtract the
        #running average cenoms and tell the user.
        for file in self.nixs_scans:
            scan_info = self.scan_info(file)
            corresponding_elastic = dir+'/'+ self.elastic_name + str.lower(scan_info[3]).split(self.nixs_name)[1]
            valid_elastic = os.path.isfile(corresponding_elastic)
            if valid_elastic:
                print("{} {}".format("Reading NIXS scan: ", file))
            else:
                print("{} {} {}".format(">> >> WARNING:", scan_info[1],"has no corresponding elastic - finding eloss by average elastic values!"))
            self.readscan_20ID(dir + '/' + file, valid_elastic)

        for file in self.wide_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading wide scan named: ", file))
            self.readscan_20ID(dir + '/' + file)

        #call function to calculate the average values over the scans - all by default
        self.average_scans(scan_numbers)
        self.get_resolutions(scan_numbers)

        #if the user asks, call function to write all info to H5 file
        if H5:
            if H5path==None:
                H5path = dir
            elif H5path:
                if os.path.isdir(H5path):
                    H5path = H5path
                else:
                    print('H5 path directory does not exist!')

            if  not sample_name:
                self.sample_name = '20ID_APS_data.H5'
                H5name = self.sample_name
            elif sample_name:
                self.sample_name = sample_name
                H5name = self.sample_name
            else:
                print('H5 sample name was not accepted')

            saveloc = H5path+'/'+H5name

            if os.path.isfile(saveloc):
                H5file = h5py.File(saveloc, "a")
            else:
                H5file = h5py.File(saveloc, "w")

            self.write_H5scanData(dir,H5file,H5name)
            print("{} {}".format("Wrote scan data to H5 file: ", saveloc))

        #let the user know the program has finished
        print('Finished Reading!')

"""To do list:
1) Check if the cenoms are close to the E0 specified in the ASCII header, and if not, do not average that scan_name
2) Make the H5 file location more interactive and allow many different samples to be read into the H5file - e.g. check if it exists and if so
write into it. DONE
3) read wide scans - Proving hard, can't do this until I've fixed the header attrbs reading, to pull out the NIXS scan region
4) Header reading and H5 attributes DONE
5) Maybe a file size check to make sure the input file isn't crazy - DONE
6) If file size is less than 5KB, then ignore it from the list - DONE
7) Read ASCII column headings to make sure that correct columns are being read DONE
8) Read ASCII headers as a saveable Info string for the USER DONE
9) WIDE SCANS!!!
10) Merge read_scans and readscan_20ID into one callable function that returns a XRStools Scan object with header attributes DONE
11) Scan_info to print the header_attrs -> scan_info to become a once stop shop for users to know the details of their scan DONE

Code to deal with the wide scans issue:
Thinking that this should be done after nixs/elastics are read in
# gives boolean index of where these conditions are True i.e. True outside of high-res region
idx = (self.scans['wide0001'].eloss > self.eloss.min())*(self.scans['wide0001'].eloss < self.eloss.max())
#this line then deletes the values in the ndarray in the high-res region
tb = np.delete(self.scans['wide0001'].eloss,np.where(idx),None)
#inserting self.eloss into the deleted gap
np.insert(tb,np.where(idx)[0][0],self.eloss)
"""
