#!/usr/bin/env python

import os
import sys
import subprocess

__all__ = ['Script_Error', 'watchrun', 'parse_cshrc', 'get_output']

class Script_Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def watchrun(cmd, err_message=None, path='.'):
#TODO: add log
    oldpath = os.getcwd()
    os.chdir(path)
    res = os.system(cmd)
    if res != 0:
        raise Script_Error(' '.join(('Error when trying to run :', \
                            cmd, 'with extra message:', str(err_message))))
    os.chdir(oldpath)

def parse_cshrc(fname):
    f = open(fname)
    for line in f:
        tokens = line.strip().split()
        if len(tokens) == 3 and tokens[0] == 'setenv':
            value = os.path.expandvars(tokens[2])
            if value.startswith('"') and value.endswith('"'):
                value = value.strip('"')
            os.environ[tokens[1]] = value

def get_output(cmd):
    return subprocess.check_output(cmd, shell=True).rstrip('\n')
    
