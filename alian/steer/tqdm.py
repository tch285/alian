from alian.steer.glob import globals

from tqdm import tqdm as tqdm_std
# import also the notebook version of tqdm
from tqdm.notebook import tqdm as tqdm_nb
tqdm_impl = tqdm_std

# determine if we are running in a notebook
try:
		from IPython import get_ipython
		if 'IPKernelApp' in get_ipython().config:
				tqdm_impl = tqdm_nb
except:
		pass

class tqdm(tqdm_impl):
	def __init__(self, *args, **kwargs):
		self.silent = kwargs.pop('silent', False)
		super(tqdm, self).__init__(*args, **kwargs)
		self.display = self.display_std
		if tqdm_impl == tqdm_nb:
				self.display = self.display_nb
		
	def display_std(self, msg=None, pos=None):
			if not globals.tqdm_silent:
					super().display(msg, pos)

	def display_nb(self, msg=None, pos=None,
									# additional signals
									close=False, bar_style=None, check_delay=True):
			if not globals.tqdm_silent:
					super().display(msg, pos, close, bar_style, check_delay)

	def write(self, s, file=None, end="\n", nolock=False):
			if not globals.tqdm_silent:
					super().write(s, file, end, nolock)
