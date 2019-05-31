#!/usr/bin/env python

'''
Plot MFM run.

Usage:
  plot [options] (<RunID> | <path>)

Options:
  -h --help    Show this screen
'''

import sys
import pickle
from docopt import docopt

from mfm import MFM
from mfm import *

def main():
    args = docopt(__doc__)

    try:
        runID = int(args['<RunID>'])
        fname = 'data/{0:03d}.mfm'.format(runID) 
    except:
        fname = args['<RunID>']
        
    if not os.path.isfile(fname):
        print("File not found: {}".format(fname))
        sys.exit()
    
    mfm = MFM(); mfm.load(fname)
    print(mfm)
    mfm.plot(PSD_seg=0.05)
    
if __name__=='__main__':
    main()
