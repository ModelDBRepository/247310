import time
import datetime
import os

class progbar(object):
    '''
    Command line progress bar.

    Parameters
    ----------
    marker : str, optional
        Progress bar marker

    Returns
    -------
    out : progbar

    Examples
    --------
    >>> from nrtlbox.utils import progbar
    >>>
    >>> pg = progbar()
    >>> for i in range(1000):
    ...     time.sleep(0.1) 
    ...     pg.update()
    ...
    '''
    def __init__(self, bar_len = 50, marker='='): 
        self.start_time = time.time()
        self.marker = marker
        self.bar_len = bar_len
        
    def update(self,pct):
        '''
        Updates and displays the progress bar

        Parameters
        ----------
        pct : float, [0,1]
            Fraction completed.
        '''

        if pct > 1: pct = 1
        if pct <= 0: pct = 1e-10
        
        elapsed = time.time() - self.start_time
        remaining = elapsed/pct - elapsed
        elapsed = str(datetime.timedelta(seconds = elapsed))
        remaining = str(datetime.timedelta(seconds = remaining))
        
        pre = '\r{:>6.2f}% {:.7s} ['.format(pct*100,elapsed)
        pre = '\r {:.7s} ['.format(elapsed)
        post = '] {:.7s} '.format(remaining)
        pct_str = ' {:>6.2f} % '.format(pct*100)

        bar = int(self.bar_len * pct) * self.marker + (self.bar_len - int(self.bar_len * pct)) * ' '

        pct_spot = int((self.bar_len - len(pct_str)) / 2)

        bar = bar[:pct_spot] + pct_str + bar[pct_spot+10:]
        
        out = pre+bar+post

        print(out,end='')

