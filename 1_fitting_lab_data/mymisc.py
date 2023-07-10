#!/usr/bin/env python

import sys, io, os
import numpy as np
import matplotlib.pylab as plt
from cycler import cycler
from contextlib import contextmanager
import ctypes
import tempfile

# some often used constants
Msun = 1.98847e33 # solar mass in g
Rsun = 6.957e10 # solar radius in cm
Lsun = 3.828e33 # solar luminosity in erg/s

G = 6.67259e-8 # gravitational constant in cgs
c = 29979245800 # cm/s

secyer = 3.1558149984e7 # seconds per year

# custom colours 
# (see https://betterfigures.org/2015/06/23/picking-a-colour-scale-for-scientific-graphics/)
col = {'black':          '#000000', 
       'orange':         '#E69F00', 
       'sky blue':       '#56B4E9', 
       'bluish green':   '#009E73', 
       'yellow':         '#F0E442', 
       'blue':           '#0072B2', 
       'vermillion':     '#D55E00', 
       'reddish purple': '#CC79A7'
      }
ls = {'-':   (), 
      '--':  (5,3), 
      ':':   (2,2), 
      '-.':  (5,3,2,3), 
      '-..': (5,3,2,3,2,3)
     }

# define my default cycler
my_colors = [col['blue'], col['vermillion'], col['sky blue'], col['orange'], col['bluish green']]
my_cycler = (cycler(color=my_colors))# + 
             #cycler(linestyle=[ls['-'], ls['--'], ls[':'], ls['-.'], ls['-..']]))

    
mpl_preamble = []
def set_plot_defaults(journal='aa', cycler=my_cycler, dpi=150, use_text_latex=True, return_midwidth=False):
    global mpl_preamble

    # page and column width of journal where plots are supposed to appear
    inches_per_pt = 1.0/72.27
    inches_per_mm = 0.03937008
    if (journal == 'aa'):
        pagewidth = 523.5307*inches_per_pt # A&A
        midwidth = 120.0*inches_per_mm
        columnwidth = 255.76535*inches_per_pt # A&A
    elif (journal == 'mnras'):
        pagewidth = 504.0*inches_per_pt # MNRAS
        midwidth = 120.0*inches_per_mm
        columnwidth = 240.0*inches_per_pt # MNRAS
        
    # enable some standard changes to plots
    fontsize = 10
    plt.rc('font', size=fontsize, family='serif')
    plt.rc('axes', prop_cycle=my_cycler, titlesize=fontsize)
    plt.rc('lines', linewidth=1, markersize=4) # a markersize of 4 roughly corresponds to using 's=13' in plt.scatter(...)
    plt.rc('xtick', top=True, direction='out')
    plt.rc('ytick', right=True, direction='out')
    plt.rc('legend', frameon=False, handletextpad=0.2)
    plt.rc('figure', figsize=(columnwidth, columnwidth * 0.6), dpi=dpi) # to enlarge inline displayed figures (saved images are vector, so this has no change)
    plt.rc('figure.constrained_layout', use=True)

    if (use_text_latex):
        # at the moment, the 'text.latex' is experimental and it can cause problems on some systems
        plt.rc('text', usetex=True)

        if (not mpl_preamble):
            mpl_preamble = plt.rcParams.get('text.latex.preamble', [])
            if (isinstance(mpl_preamble, str)):
                mpl_preamble = [mpl_preamble]
            mpl_preamble = mpl_preamble +   [
                                                r'\usepackage{amsmath}', 
                                                r'\usepackage{txfonts}', 
                                                r'\newcommand{\msun}{\ensuremath{\mathrm{M}_{\odot}}}', 
                                                r'\newcommand{\lsun}{\ensuremath{\mathrm{L}_{\odot}}}', 
                                                r'\newcommand{\rsun}{\ensuremath{\mathrm{R}_{\odot}}}', 
                                                r'\newcommand{\zsun}{\ensuremath{\mathrm{Z}_{\odot}}}'
                                            ]

        plt.rc('text.latex', preamble=mpl_preamble)
                                    
    if (return_midwidth):
        return pagewidth, columnwidth, midwidth
    else:
        return pagewidth, columnwidth
    

def disable_stdout(f):
	# taken from Arepo Snaputils: stellar_ics/tools.py
	def fun(*arg, **kwarg):
		stdout = sys.stdout
		sys.stdout = io.StringIO()
		result = None
		try:
			result = f(*arg, **kwarg)
		finally:
			sys.stdout = stdout
		return result
	return fun


# context manager solution to also disable stdout from C libraries
# taken from https://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/
#
# Usage example:
# f = io.BytesIO()
# 
# with stdout_redirector(f):
#     print('foobar')
#     print(12)
#     libc.puts(b'this comes from C')
#     os.system('echo and this is from echo')
# print('Got stdout: "{0}"'.format(f.getvalue().decode('utf-8')))
#
libc = ctypes.CDLL(None)
c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')

@contextmanager
def stdout_redirector(stream):
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stdout_fd)


def execute_notebook(nbfile):
    import nbformat

    with (io.open(nbfile, encoding='utf8')) as f:
        nb = nbformat.read(f, as_version=4)

    ip = get_ipython()

    for cell in nb.cells:
        if (cell.cell_type != 'code'):
            continue
        ip.run_cell(cell.source)

