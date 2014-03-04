"""
To use this:

    sudo pip install oct2py
"""

import os.path

import oct2py as op

def solvedbo(u_q):
    oct_session = op.Oct2Py()
    oct_session.addpath(os.path.join(os.path.expanduser('~'), 'work/github/phenology-two-trait-migratory-bird'))

    oct_session.run('[x_cV, yzV, nV, z_n] = solvedbowrapper(%f)' % u_q)

    return { 'x_cV': oct_session.get('x_cV'),
             'yzV':  oct_session.get('yzV'),
             'nV':   oct_session.get('nV'),
             'z_n':  oct_session.get('z_n'),
           }
