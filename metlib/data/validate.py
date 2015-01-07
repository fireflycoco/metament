#!/usr/bin/env python
import numpy as np

def rec_validate(rec, by_fields=None):
    fields_to_check = []
    for f, (dt, offset) in rec.dtype.fields.iteritems():
        if np.issubdtype(dt, np.float):
            print f, dt
