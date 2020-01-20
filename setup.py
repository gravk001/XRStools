#!/usr/bin/env python
#-*- coding: utf-8 -*-
#
#    Project: xrstools
#

"""
Installer script for xrstools
"""

__authors__   = ["Ch.J. Sahle"]
__contact__   = "christoph.sahle@helsinki.fi"
__license__   = "a voir"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France "
__date__      = "2014-06-03"
__status__    = "beta"


import os, sys, glob, shutil, ConfigParser, platform
from distutils.core import setup, Extension, Command
from numpy.distutils.misc_util import get_numpy_include_dirs
from distutils.sysconfig import get_python_lib
import string
import stat


try:
    import numpy
except ImportError:
    text = "You must have numpy installed.\n"
    raise Exception(text)

try:
    from Cython.Distutils import build_ext
    CYTHON = True
except ImportError:
    CYTHON = False

if CYTHON:
    cython_c_ext = ".pyx"
else:
    cython_c_ext = ".cpp"
    from distutils.command.build_ext import build_ext

data_files =  glob.glob("data/*.dat")+glob.glob("data/chitables/*.dat")

from distutils.command.install_data import install_data

class smart_install_data(install_data):
    def run(self):
        global  version
        install_cmd = self.get_finalized_command('install')
        self.install_dir = os.path.join(getattr(install_cmd, 'install_lib'),"xrstools"+version,"..","..","..","..","share","xrstools","data")
        print("DATA to be installed in %s" % self.install_dir)

        install_lib_dir = os.path.join(getattr(install_cmd, 'install_lib'),"xrstools"+version )
        if not os.path.exists(self.install_dir):
            os.makedirs(self.install_dir)
        print " QUI "
        for filein in glob.glob('xrstools/*ui'):
            filedest = os.path.join(install_lib_dir, os.path.basename(filein))
            print "cp %s %s"%(filein,filedest)
            os.system("cp %s %s"%(filein,filedest))

        return install_data.run(self)

###########################################################################################################
from distutils.command.install_scripts import install_scripts
class smart_install_scripts(install_scripts):
    def run (self):
        global  version

        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_scripts')
        self.install_lib_dir = getattr(install_cmd, 'install_lib')

        self.outfiles = []
        for filein in glob.glob('xrstools/scripts/*'):
            if filein in glob.glob('xrstools/scripts/*~'): continue

            print "INSTALLO ", filein

            if not os.path.exists(self.install_dir):
                os.makedirs(self.install_dir)

            filedest = os.path.join(self.install_dir, os.path.basename(filein+version))
            print os.path.basename(filein+version)

            if os.path.exists(filedest):
                os.remove(filedest)

            text = open(filein, 'r').read()
            text=string.replace(text,"import xrstools", "import xrstools"+version)
            text=string.replace(text,"python", sys.executable)

            f=open(filedest, 'w')
            f.write(text)
            f.close()
            self.outfiles.append(filedest)

        if os.name == 'posix':
            # Set the executable bits (owner, group, and world) on
            # all the scripts we just installed.
            for file in self.get_outputs():
                if self.dry_run:
                    pass
                else:
                    mode = ((os.stat(file)[stat.ST_MODE]) | 0555) & 07777
                    os.chmod(file, mode)

if sys.platform == "win32":
	pass
    ## This is for mingw32/gomp?
    #data_files[0][1].append(os.path.join("dll", "pthreadGC2.dll"))
    #root = os.path.dirname(os.path.abspath(__file__))
    #tocopy_files = []
    #script_files = []
    #for i in os.listdir(os.path.join(root, "scripts")):
    #    if os.path.isfile(os.path.join(root, "scripts", i)):
    #        if i.endswith(".py"):
    #            script_files.append(os.path.join("scripts", i))
    #        else:
    #            tocopy_files.append(os.path.join("scripts", i))
    #for i in tocopy_files:
    #    filein = os.path.join(root, i)
    #    if (filein + ".py") not in script_files:
    #        shutil.copyfile(filein, filein + ".py")
    #        script_files.append(filein + ".py")
else:
    script_files = glob.glob("./xrstools/scripts/*")

version = '' #[eval(l.split("=")[1]) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "xrstools", "__init__.py")) if l.strip().startswith("version")][0]
# version = '0.0' #[eval(l.split("=")[1]) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "xrstools", "__init__.py")) if l.strip().startswith("version")][0]

# We subclass the build_ext class in order to handle compiler flags
# for openmp and opencl etc in a cross platform way
translator = {
        # Compiler
            # name, compileflag, linkflag
        'msvc' : {
            'openmp' : ('/openmp', ' '),
            'debug'  : ('/Zi', ' '),
            'OpenCL' : 'OpenCL',
            },
        'mingw32':{
            'openmp' : ('-fopenmp', '-fopenmp'),
            'debug'  : ('-g', '-g'),
            'stdc++' : 'stdc++',
            'OpenCL' : 'OpenCL'
            },
        'default':{
            'openmp' : ('-fopenmp', '-fopenmp'),
            'debug'  : ('-g', '-g'),
            'stdc++' : 'stdc++',
            'OpenCL' : 'OpenCL'
            }
        }

cmdclass = {}

class build_ext_tds2el(build_ext):
    ## this is for templating
    def build_extensions(self):

        if self.compiler.compiler_type in translator:
            trans = translator[self.compiler.compiler_type]
        else:
            trans = translator['default']



        for e in self.extensions:
            e.extra_compile_args = [ trans[a][0] if a in trans else a
                                    for a in e.extra_compile_args]
            e.extra_link_args = [ trans[a][1] if a in trans else a
                                 for a in e.extra_link_args]

            if    e.libraries !=[ "fftw3f"] and e.libraries !=[ "mpi_cxx"]:
                e.libraries = filter(None, [ trans[a] if a in trans else None
                                             for a in e.libraries])



            # If you are confused look here:
            # print e, e.libraries
            # print e.extra_compile_args
            # print e.extra_link_args
        build_ext.build_extensions(self)

# cmdclass['build_ext'] = build_ext_tds2el



class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        import sys, subprocess
        os.chdir("test")
        errno = subprocess.call([sys.executable, 'test_all.py'])
        if errno != 0:
            raise SystemExit(errno)
        else:
            os.chdir("..")
cmdclass['test'] = PyTest

#######################
# build_doc commandes #
#######################

try:
    import sphinx
    import sphinx.util.console
    sphinx.util.console.color_terminal = lambda: False
    from sphinx.setup_command import BuildDoc
except ImportError:
    sphinx = None

if sphinx:
    class build_doc(BuildDoc):

        def run(self):

            # make sure the python path is pointing to the newly built
            # code so that the documentation is built on this and not a
            # previously installed version

            build = self.get_finalized_command('build')
            sys.path.insert(0, os.path.abspath(build.build_lib))

            # Build the Users Guide in HTML and TeX format
            for builder in ('html', 'latex'):
                self.builder = builder
                self.builder_target_dir = os.path.join(self.build_dir, builder)
                self.mkpath(self.builder_target_dir)
                builder_index = 'index_{0}.txt'.format(builder)
                BuildDoc.run(self)
            sys.path.pop(0)
    cmdclass['build_doc'] = build_doc

ext_modules = []

#############################################################################

jn = os.sep.join

### TEMPLATE#################
c_sorgenti_tmp = ["spotpicker_cy"+cython_c_ext, "spotpicker.cc"]
c_sorgenti = [jn(['tds2el', 'tds2el_c', tok]) for tok in c_sorgenti_tmp]
depends_tmp = [ "spotpicker.h"]
depends = [jn(['tds2el', 'tds2el_c', tok]) for tok in depends_tmp]
cython_ext2 = Extension(name="spotpicker_cy",
                        sources=c_sorgenti,
                        depends=depends,
                        libraries=[ "mpi_cxx"],
                        include_dirs=get_numpy_include_dirs()+["/usr/include/openmpi/"],
                        language="c++",             # generate C++ code
                        extra_compile_args=['-fopenmp'] ,
                        extra_link_args=['-fopenmp'] ,
                        )
## ext_modules.append(cython_ext2)
##########################################

cmdclass['install_scripts'] = smart_install_scripts
cmdclass['install_data'] = smart_install_data

# setup(name          = 'xrstools',
#       version       = version,
#       author        = "Christoph J. Sahle",
#       author_email  = "christoph.sahle@esrf.fr",
#       description   = 'An XRS data analysis package.',
#       url           = "christoph.sahle@esrf.fr",
#       download_url  = "",
#       ext_package   = "xrstools"+version,
#       scripts       = script_files,
#       packages      = ["xrstools"+version],
#       package_dir   = {"xrstools"+version:"xrstools" },
#       package_data  = {"xrstools"+version: ['data/*.dat','data/chitables/*.dat','data/refl_database/*.dat','docs/*.txt','examples/*.*']},
#       test_suite    = "test",
#       cmdclass      = cmdclass,
#       ext_modules   = ext_modules
# )






setup(name='xrstools',
      version=version,
      author="Christoph Sahal",
      author_email="",
      description='descrizione da farsi qui',
      url="christoph.sahle@esrf.fr",
      download_url="",
      ext_package="xrstools"+version,
      scripts=glob.glob("./xrstools/scripts/*"),
      # ext_modules=[Extension(**dico) for dico in ext_modules],
      packages=["xrstools"+version],
      package_dir={"xrstools"+version:"xrstools" },
      test_suite="test",
      cmdclass=cmdclass,
      data_files=data_files,
      ext_modules=ext_modules
      )

# try:
#     import pyopencl
# except ImportError:
#     print("""sprsaocl can use pyopencl to run on parallel accelerators like GPU; this is an optional dependency.
# This python module can be found on:
# http://pypi.python.org/pypi/pyopencl
# """)
