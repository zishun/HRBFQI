import sys
sys.path.append('./build/')

import numpy as np
from pyHRBFQI import hrbfqi

def load_pwn(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    n = int(lines[0])
    pnts = np.zeros((n,3))
    nmls = np.zeros((n,3))
    for i in range(n):
        pnts[i,:] = [float(x) for x in lines[i+1].split()]
        nmls[i,:] = [float(x) for x in lines[n+i+1].split()]
    return pnts, nmls

def write_obj(fn, verts, faces):
    with open(fn, 'w') as f:
        for v in verts:
            f.write('v %.8f %.8f %.8f\n'%(v[0],v[1],v[2]))
        for ff in faces+1:
            f.write('f %d %d %d\n'%(ff[0],ff[1],ff[2]))

pnts, nmls = load_pwn('../Bin/Data/torus.pwn')
mesh_v, mesh_f = hrbfqi(pnts, -nmls, False, False, 1.2, 50, 0.1, 0, 0, True)
write_obj('./build/output.obj', mesh_v, mesh_f)
