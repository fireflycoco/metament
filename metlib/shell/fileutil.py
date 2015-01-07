#!/usr/bin/env python
import os
from glob import glob
import re
import shutil
from metlib.shell.script_helper import get_output

__all__ = ['filesize', 'sorted_walk', 'list_all_file',  
            'force_rm', 'force_makedirs', 'expand_path', 'CB', 'get_rel_path', 'find_link_orig', 
            'strip_ext', 'sub_ext', 'get_ext', 
            'LS', 'LS_R', 'CD', 'P', 'DIRNAME', 'BASENAME', 
            'RM', 'CP', 'MKDIR', 'MV']

def filesize(f):
    """Return the size of f in bytes"""
    if isinstance(f, file):
        now_pos = f.tell()
        f.seek(0, 2)
        size = f.tell()
        f.seek(now_pos)
        return size
    else:
        return os.stat(f).st_size

def sorted_walk(top, **kwarg):
    """returns the sorted result in a list consisting os.walk's 3-tuple: (dirpath, dirnames, filename). 
    kwargs: os.path's kwargs, i.e. topdown=True[, onerror=None[, followlinks=False]]"""
    if 'followlinks' not in kwarg:
        kwarg['followlinks'] = True
    w = os.walk(top, **kwarg)
    res = list(w)
    res.sort()
    for stuff in res:
        stuff[1].sort()
        stuff[2].sort()
    return res

def list_all_file(top='.', fname_pattern=r'.*', dir_pattern=r'.*', \
                    ignore_hidden=True, recursive=True, **kwarg):
    """returns a list of filenames that matches the 2 regex patterns, 
    kwargs: os.path's kwargs, i.e. topdown=True[, onerror=None[, followlinks=False]]"""
    try:
        newtop = os.path.expanduser(os.path.expandvars(top))
    except Exception:
        newtop = top
    if 'followlinks' not in kwarg:
        kwarg['followlinks'] = True
    if recursive:
        walktuplelist = sorted_walk(newtop, **kwarg)
    else:
        fnames = sorted(os.listdir(newtop))
        walktuplelist = [('', [], fnames)]
    filelist = []
    for walktuple in walktuplelist:
        d = walktuple[0]
        if re.search(dir_pattern, d):
            if ignore_hidden:
                sorted_dirs = sorted(d.split(os.path.sep))
                sorted_dirs = [subd for subd in sorted_dirs if subd not in ('', '.', '..')]
                if len(sorted_dirs) > 0 and sorted_dirs[0].startswith('.'):
                    continue
            for fname in walktuple[2]:
                if re.search(fname_pattern, fname):
                    if ignore_hidden:
                        if fname.startswith('.') or fname.endswith('~'):
                            continue
                    filelist.append(os.path.join(d, fname))
    return filelist

def LS(top='.', **kwarg):
    """wrapping for list_all_file() with recursive default as False.
    kwargs are the same as list_all_file, i.e.,
        fname_pattern, dir_pattern, ignore_hidden, recursive, followlinks, etc.
    """
    recursive = kwarg.pop('recursive', False)
    return list_all_file(top, recursive=recursive, **kwarg)
    
LS_R = list_all_file

def CD(path=None):
    """CD changes dir
    """
    if path is None:
        path = P('~')
    else:
        path = P(path)
    return os.chdir(path)

def expand_path(path):
    """expand ~ and $ in the path"""
    return os.path.expanduser(os.path.expandvars(path))

def CB(tag):
    """Using the metamet/shell/cb (short for clipboard) to get recorded text snippets, e.g., long file names, etc."""
    return get_output('cb %s' % tag)

def force_rm(fname, regex=False):
    """force to rm fname, no matter whether fname is file, dir or link"""
    # TODO, add regex support. and make it more robust.
    fnames = glob(fname)
    for fn in fnames:
        try:
            if os.path.islink(fn):
                os.unlink(fn)
            elif os.path.isdir(fn):
                shutil.rmtree(fn)
            else:
                os.remove(fn)
        except Exception as e:
            print e

def force_makedirs(dirname, rm_exist_dir=False):
    """force to make dir, no matter whether it exists"""
    #TODO: make it more robust
    if dirname in ['', '.', '..', './', '../', '/']:
        return
    orig = dirname
    if os.path.islink(orig):
        orig = find_link_orig(dirname)
    if os.path.isfile(orig):
        if not rm_exist_dir:
            raise RuntimeError('Cannot makedirs: %s is file' % dirname)
    if rm_exist_dir:
        try:
            force_rm(dirname)
        except Exception as e:
            pass
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != 17:
            raise e

def get_rel_path(path, base):
    """get relative path, e.g., get_rel_path('abc/de/fg', 'abc') => 'de/fg'
    """
    lb = len(base)
    assert path[:lb] == base
    if len(path) == lb:
        rel_path = ''
    elif path[lb] == '/':
        rel_path = path[lb+1:]
    else:
        rel_path = path[lb:]
    return rel_path

def find_link_orig(path, max_depth=99):
    """Try to find the orig of a link."""
    count = 0
    while count < max_depth:
        if os.path.islink(path):
            path = os.readlink(path)
        else:
            return path
        count += 1
    return path

P = expand_path
RM = force_rm
MKDIR = force_makedirs
DIRNAME = os.path.dirname
BASENAME = os.path.basename
MV = shutil.move

def CP(src, dest, symlinks=True, ignore=None):
    true_dest = find_link_orig(dest)
    if os.path.isdir(src):
        if os.path.isdir(true_dest):
            dest = os.path.join(dest, os.path.basename(src))
        else:
            RM(dest)
        return shutil.copytree(src, dest)
    else:
        if not os.path.isdir(true_dest):
            RM(dest)
        return shutil.copy(src, dest)

def strip_ext(path):
    """strips .ext from a path"""
    return os.path.splitext(path)[0]

def get_ext(path):
    """get .ext from a path"""
    return os.path.splitext(path)[1]

def sub_ext(orig, new_ext):
    """sub .ext with a new one"""
    if not new_ext.startswith('.'):
        new_ext = '.' + new_ext
    return strip_ext(orig) + new_ext
