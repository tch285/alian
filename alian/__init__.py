import sys
import os

from yasp import GenericObject

class AliAnSettings(GenericObject):
	def __init__(self, **kwargs):
		super(AliAnSettings, self).__init__(**kwargs)


alian_settings = AliAnSettings(src_path=os.path.dirname(os.path.abspath(__file__)), path=os.environ['ALIAN_DIR'])

# here load all the ${ALIAN_DIR}/include/*.hh files with cppyy and the libraries found under ..${ALIAN_DIR}/lib
# use cppyyhelper to do this

from yasp import find_files

headers = find_files(alian_settings.path + '/include', '*.hh')
packs = ['alian']
libs = set(['alian_' + os.path.basename(os.path.dirname(h)) for h in headers])

from yasp.cppyyhelper import YaspCppyyHelper
YaspCppyyHelper().load(packs, libs, headers)
print(YaspCppyyHelper())

# import heppyy
# alian = heppyy.load_cppyy('alian')
