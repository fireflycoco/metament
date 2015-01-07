#!/usr/bin/env python

# misc.py

import os, sys
import re

#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
from numpy import array
import collections
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset
from metlib.misc.datatype import Null
from metlib.misc.datatype import isseq

__all__ = ['struni', 'grep', 'str2list', 'get_sys_argv', 'parse_bool', 'Setter']

def struni(obj):
    """ return str(obj) if possible, else return unicode(obj).
    """
    try:
        res = str(obj)
        return res
    except UnicodeEncodeError as e:
        res = unicode(obj)
        return res

def grep(pattern, seq, flags=0):
    """grep greps patterns from seqs.
Parameters:
    pattern: regex pattern str or compiled regex pattern or func. 
    seq: seq of str or any other stuffs, etc. If pattern is not func and seq[0] is not str/unicode, seq will be converted using str() first.
    flags: 'ILMSUX' or number of re.I|re.L|re.M|re.S|re.U|re.X .
        I: ignore case;
        L: locale dependent;
        M: multi-line;
        S: dot matches all;
        U: unicode dependent;
        X: verbose.
Returns:
    A filtered list.
    """
    if isinstance(flags, (str, unicode)):
        true_flags = 0
        flags = flags.upper()
        for letter in flags:
            if letter in 'ILMSUX':
                true_flags += re.__dict__[letter]
    else:
        true_flags = flags
    if isinstance(pattern, (str, unicode, type(re.compile(r'')) ) ):
        res = []
        for item in seq:
            if re.search(pattern, struni(item), true_flags):
                res.append(item)
    elif callable(pattern):
        res = filter(pattern, seq)

    return res


def str2list(s, pattern=',|;|:|#|\||\s+'):
    return re.split(pattern, s)

def get_sys_argv(argnames=[], optional_argnames=[], keyword_argnames=[], desc=''):
    """get_sys_argv parse sys.argv according to a list of argnames. 
    Parsed args goes directly into globals().
    If not succeed, print a usage prompt and exit.

Parameters
----------
argnames: a list of arg names or tuples of (arg name, convert function).
optional_argnames: like argnames, but optional. It's best to specify default values for these optional args before calling get_sys_argv(), as shown in the example; otherwise the default value will be Null.
keyword_argnames: optional keyword args (e.g., kw1=..., kw2=...). If not specified, *NO* variable will be generated. So it's best to specify non-optional keyword values before calling get_sys_argv(). 
desc: additional descriptions.

Return
------
a list of remaining args.
And two extra dicts will be generated: _sys_argv_args, _sys_argv_kwargs

Example
-------
>>> job = 'sinner'
>>> status = 'not redeemed yet'
>>> get_sys_args( ['name', ('age', int)], ['job', 'status'])

    """
    class ExceptionHelp(Exception):
        pass

    nonkw_argv = []
    kw_argv = []
    for a in sys.argv[1:]:
        if '=' in a:
            kw_argv.append(a)
        else:
            nonkw_argv.append(a)
    kwargconvd = dict(keyword_argnames)
    argdict = dict()
    kwargdict = dict()
    true_argnames = []
    argconvd = dict()
    n_must = len(argnames)
    for a in argnames + optional_argnames:
        if isinstance(a, (str, unicode)):
            true_argnames.append(str(a))
        elif isinstance(a, (tuple, list)) \
                and len(a) == 2 \
                and isinstance(a[0], (str, unicode)) \
                and callable(a[1]):
            true_argnames.append(str(a[0]))
            argconvd[a[0]] = a[1]
    try:
        if len(sys.argv) == 2 and sys.argv[1] == "--help":
            raise ExceptionHelp()
        for i, argname in enumerate(true_argnames):
            if i >= len(nonkw_argv) and i >= n_must:
                for argn in true_argnames[i:]:
                    if argn not in sys.modules['__main__'].__dict__:
                        sys.modules['__main__'].__dict__[argn] = Null
                break
            value = nonkw_argv[i]
            if argname in argconvd:
                value = argconvd[argname](value)
            argdict[argname] = value
        for argstr in kw_argv:
            argname, value = argstr.split('=', 1)
            if argname in kwargconvd:
                value = kwargconvd[argname](value)
            kwargdict[argname] = value
        sys.modules['__main__'].__dict__.update(argdict)
        sys.modules['__main__'].__dict__.update(kwargdict)
        sys.modules['__main__'].__dict__['_sys_argv_args'] = argdict
        sys.modules['__main__'].__dict__['_sys_argv_kwargs'] = kwargdict
    except Exception as e:
        print "Usage:"
        print "    %s" % sys.argv[0], 
        for argn in true_argnames[:n_must]:
            print argn,
        if len(optional_argnames) > 0:
            print '[',
            for argn in true_argnames[n_must:]:
                print argn,
            print ']',
        if len(kwargconvd) > 0:
            print '{',
            for kwargn in kwargconvd:
                print '%s=' % kwargn,
            print '}',
        print
        for argn in true_argnames:
            if argn in argconvd:
                print "        %s will be converted by %s" % (argn, argconvd[argn])
        for kwargn, convfunc in kwargconvd.iteritems():
            if convfunc not in (str, unicode):
                print "        %s will be converted by %s" % (kwargn, convfunc)
        print desc
        if not isinstance(e, ExceptionHelp):
            print "Error:"
            print "    %s:" % argname, e
        sys.exit(1)
    other_args = nonkw_argv[len(true_argnames):]

    return other_args

def parse_bool(s):
    def _tobool(ss):
        ss = ss.upper()
        if ss in ['FALSE', 'F', 'NO', 'N', '0', '']:
            return False
        elif ss in ['TRUE', 'T', 'YES', 'Y', '1']:
            return True
        else:
            return bool(ss)

    if isinstance(s, (str, unicode)):
        if s.upper() in ['TRUE', 'FALSE', 'T', 'F', 'YES', 'NO', 'Y', 'N', '1', '0', '']:
            return _tobool(s)
        elif set(s.upper()) <= set(['T', 'F', '0', '1', 'Y', 'N']):
            return np.array([_tobool(ss) for ss in s])
        else: 
            return bool(s)
    elif isinstance(s, (collections.Sequence, np.ndarray)):
        return np.array([_tobool(ss) for ss in s])
    else:
        return bool(s)

class Setter(object):
    """Setter is a base-class for loading and saving settings.
    Example
    -------
    >>> class HasSetter(Setter):
    ...     def __init__(self, a, b):
    ...         self.a = a
    ...         self.b = b
    ... o = HasSetter(3, 5)
    ... o.save_setting('settings.rc', prefix='set_o.', attrs=['a', 'b'])

    settings.rc will be like:
        set_o.a = 3
        set_o.b = 5

    >>> o2 = HasSetter(0, 0)
    ... print o2.a, o2.b
    <<< 0 0
    
    >>> o2.load_setting('settings.rc', prefix='set_o.')
    ... print o2.a, o2.b
    <<< 3 5

    """ 
    def load_setting(self, setting, prefix=''):
        if isseq(setting):
            setting_list = setting
        elif isinstance(setting, (str, unicode)):
            setting_file = open(setting)
            setting_list = []
            current_lines = []
            for l in setting_file:
                l = l.rstrip()
                if not l:
                    continue
                current_lines.append(l)
                if not l.endswith(','):
                    joined = '\n'.join(current_lines)
                    setting_list.append(joined)
                    current_lines = []
            if current_lines:
                setting_list.append('\n'.join(current_lines))
        for setting_item in setting_list:
            if setting_item.startswith(prefix):
                try:
                    setting_item = setting_item.lstrip(prefix)
                    cmd = 'self.' + setting_item
                    exec cmd
                except Exception as e:
                    pass

    def save_setting(self, setting_file, filemode='w', prefix='', attrs=None):
        setting_list = self.get_setting(prefix, attrs)
        with open(setting_file, filemode) as outfile:
            for setting_item in setting_list:
                outfile.write(setting_item)
                outfile.write('\n')

    def get_setting(self, prefix='', attrs=None):
        res = []
        if attrs is None:
            for attr, value in self.__dict__.iteritems():
                res.append('%s%s = %s' % (prefix, attr, repr(value)))
        else:
            for attr in attrs:
                res.append('%s%s = %s' % (prefix, attr, repr(eval('self.%s' % attr))))
        return res

if __name__ == '__main__':
    l = ['abcde', 'asldkfj', 'sdjfowij', 'sdfoijw', '1243450', '1204023']
    print grep(re.compile('sd'), l)
    print grep(r"^\d*$", l)
    print grep(lambda s: s.startswith('a'), l)
  

