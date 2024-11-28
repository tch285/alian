import sys

class ColorS(object):
	def str(*args):
		_s = ' '.join([str(s) for s in args])
		return _s
	def red(*s): 
		return '\033[91m{}\033[00m'.format(ColorS.str(*s))
	def green(*s): 
		return '\033[92m{}\033[00m'.format(ColorS.str(*s))
	def yellow(*s): 
		return '\033[93m{}\033[00m'.format(ColorS.str(*s))
	def blue(*s): 
		return '\033[34m{}\033[00m'.format(ColorS.str(*s))
	def light_purple(*s): 
		return '\033[94m{}\033[00m'.format(ColorS.str(*s))
	def purple(*s): 
		return '\033[95m{}\033[00m'.format(ColorS.str(*s))
	def cyan(*s): 
		return '\033[96m{}\033[00m'.format(ColorS.str(*s))
	def light_gray(*s): 
		return '\033[97m{}\033[00m'.format(ColorS.str(*s))
	def no_color(*s): 
		return '\033[00m{}\033[00m'.format(ColorS.str(*s))
	def black(*s):
		return '\033[98m{}\033[00m'.format(ColorS.str(*s))
	def __init__(self):
		pass
	# credit: https://www.geeksforgeeks.org/print-colors-python-terminal/
	# Python program to print 
	# colored text and background 
	def print_format_table():
		""" 
		prints table of formatted text format options 
		"""
		for style in range(8): 
			for fg in range(30, 38): 
				s1 = '' 
				for bg in range(40, 48): 
					format = ';'.join([str(style), str(fg), str(bg)]) 
					s1 += '\x1b[%sm %s \x1b[0m' % (format, format) 
				print(s1) 
			print('\n') 

def pwarning(*args, file=sys.stderr):
	print(ColorS.yellow('[w]', *args), file=file)

def pdebug(*args, file=sys.stderr):
	print(ColorS.purple('[d]', *args), file=file)

def perror(*args, file=sys.stderr):
	print(ColorS.red('[e]', *args), file=file)

def pinfo(*args, file=sys.stdout):
	print(ColorS.green('[i]', *args), file=file)

def pindent(*args, file=sys.stdout):
	print(ColorS.no_color('   ', *args), file=file)

class CursorSpin(object):
	_cpos = 0
	_cursor = '\\|/-'
	def __init__(self):
		sys.stdout.write(' {}\r'.format(CursorSpin._cursor[CursorSpin._cpos]))
		sys.stdout.flush()
		CursorSpin._cpos = CursorSpin._cpos + 1
		if CursorSpin._cpos >= len(CursorSpin._cursor):
			CursorSpin._cpos = 0
