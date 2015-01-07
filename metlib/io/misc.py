import pickle

__all__ = ['savepickle', 'loadpickle']

def savepickle(fname, obj):
    outf = open(fname, 'w')
    pickle.dump(obj, outf)
    outf.close()

def loadpickle(fname):
    infile = open(fname)
    obj = pickle.load(infile)
    infile.close()
    return obj