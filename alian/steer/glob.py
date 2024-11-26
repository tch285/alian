from heppyy import GenericObject
# singleton class for storing unique analysis names

class Globals(GenericObject):

		__instance = None

		def __new__(cls):
				if Globals.__instance is None:
						Globals.__instance = object.__new__(cls)
						Globals.__instance.names = []
				return Globals.__instance
  
globals = Globals()
globals.tqdm_silent = False
globals.debug = False
globals.verbose = False
globals.log = None
