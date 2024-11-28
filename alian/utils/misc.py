def is_iterable(o):
	result = False
	try:
		tmp_iterator = iter(o)
		result = True
	except TypeError as te:
		result = False
	return result

def is_iterable_not_string(o):
	result = False
	try:
		tmp_iterator = iter(o)
		result = True
	except TypeError as te:
		result = False
	if type(o) == str:
		result = False
	return result

# think about thread safe implementation
# use unique file names... for example?
class UniqueString(object):
	locked_strings = []
	def __init__(self, base=None):
		self.base = base

	def _unique(base=None):
		i = 0
		retstring = base
		if retstring is None:
			retstring = 'UniqueString_0'
		else:
			retstring = '{}_{}'.format(str(base), i)
		while retstring in UniqueString.locked_strings:
			retstring = '{}_{}'.format(retstring.split('_')[0], i)
			i = i + 1
		UniqueString.locked_strings.append(retstring)
		return retstring

	def str(self, base=None):
		if base:
			self.base = base
		return UniqueString._unique(self.base)

	def str(base=None):
		return UniqueString._unique(base)


class NoneSetWrapper(object):
	def __init__(self, name):
		self.name = name
	def __getattr__(self, key):
		try:
			return self.__dict__[key]
		except:
			print(ColorS.red('[w] {} : {} attribute is not known'.format(self.name, key)), file=sys.stderr)
			self.__setattr__(key, None)
		return None
	def description(self):
		return 'NoneSetWrapper named {}'.format(self.name)

class NoneSetWrappers(object):
	_instance = None
	def __init__(self):
		self.wrappers = {}
	def get(self, name):
		try:
			return self.wrappers[name]
		except KeyError:
			self.wrappers[name] = NoneSetWrapper(name)
		return self.wrappers[name]
	def instance():
		if NoneSetWrappers._instance is None:
			NoneSetWrappers._instance = NoneSetWrappers()
		return NoneSetWrappers._instance

class Type(object):
	_float = type(0.)
	_int = type(0)
	_list = type([0,1])
	_tuple = type((0,1))
	def __init__(self):
		pass
	def is_float(x):
		return (float == type(x))
	def is_int(x):
		return (int == type(x))
	def is_list(x):
		return (list == type(x))
	def is_tuple(x):
		return (tuple == type(x))
	def is_dict(x):
		return (dict == type(x))
