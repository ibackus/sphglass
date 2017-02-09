# -*- coding: utf-8 -*-
"""
Contains functions for generating an SPH glass in a periodic box.  To generate
a glass, see glassBox.

This package requires diskpy, ChaNGa, and pynbody.


Created on Wed Mar 16 17:37:10 2016

@author: ibackus
"""
import shutil
import os
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import diskpy
from diskpy.ICgen.ICgen_utils import changa_command, changa_run
import itertools

# Constants
defaultparam = 'glassdefaults.param'
directory = os.path.split(os.path.abspath(__file__))[0]
defaultparam = os.path.join(directory, defaultparam)
kB = SimArray(1.0, 'k')

# Set up default params if required (ie the first time this package is run)
if not os.path.exists(defaultparam):
    
    defaults = os.path.join(directory, '.defaultparam')
    shutil.copyfile(defaults, defaultparam)
    print 'Setting up default params...saved to ' + defaultparam
    
def _loadDefaults():
    """
    Load default .param file
    """
    return diskpy.utils.configparser(defaultparam, 'param')
    
def filenames():
    """
    Return default filenames
    """
    param = _loadDefaults()
    inFile = param['achInFile']
    outPrefix = param['achOutName']
    
    return inFile, outPrefix

def glassBox(n, shape=[1,1,1], changaPreset='default', verbose=False,
             fulloutput=False, usegrid=False, randomness=1., nSmooth=32, 
             runTimeScale=1., dampingScale=1., extraPars={}):
    """
    Generates an sph glass in a box with periodic boundary conditions using 
    ChaNGa.  The procedure is:
        1. Generate random particle positions in a box
        2. Create a tipsy snapshot with only gas particles
        3. Time evolve in a periodic box with no gravity and lots of artificial
            viscosity.
    
    ND-SPH is supported for SPH in dimensions 1, 2, and 3 IF ChaNGa has been
    properly compiled with NDPSH.  The number of dimensions to run in is 
    from the length of shape.
    
    Parameters
    ----------
    n : int or list/array-like
        Number of particles (if int) or grid resolution [nx, ny, nz] along each
        axis (only if usegrid=True) 
    shape : array-like
        Shape of the box (x, y, z).  The box will be centered around the origin.
    changaPreset : str
        ChaNGa preset to use (see diskpy for info on the ChaNGa presets)
    verbose : bool
        Verbosity, true or false
    fulloutput : bool
        If True, all the snapshots for each time step during the time evolution
        will be output.
        This is useful especially with long boxes where waves can form.
    usegrid : bool
        Use a grid to seed intial positions.  The particles will be randomly 
        shifted around the grid locations
    randomness : float
        If usegrid=True, specifies by what fraction of the grid spacing
        particles will be randomly shifted
    nSmooth : int
        Number of neighbors to use for SPH
    runTimeScale : float
        Factor to increase ChaNGa run time by.  If you think your glasses are
        not getting fully settled into a glass state, try increasing this
        number.
    dampingScale : float
        Factor to increase the damping force in ChaNGa.
    extraPars : dict
        param dict defining params to override the default ChaNGa runtime
        params defined here.
    
    
    Returns
    -------
    f : SimSnap (see pynbody)
        Snapshot containing the generated particle positions.  The glass
        positions can be accessed using f['pos'].  Density is also available for
        all the particles in f['rho'].  The target density is 1
    
    Notes
    -----
    The snapshot will saved to glass.std
    
    If you find the box is not sufficiently glassy, you can time evolve it
    again by running reglassify() which will run glass.std again.
    """
    # Generate snapshot with random positions
    snap = boxSnap(n, shape, usegrid, randomness)
    param = makeParam(snap, shape, fulloutput, runTimeScale, dampingScale)
    # Save snapshot and param
    paramname = param['achOutName'] + '.param'
    ICname = param['achInFile']
    param['nSmooth'] = nSmooth
    param.update(extraPars)
    diskpy.utils.configsave(param, paramname)
    snap.write(fmt=pynbody.tipsy.TipsySnap, filename=ICname)
    # Run ChaNGa to make a glass
    f = runchanga(paramname, changaPreset, verbose, fulloutput)
    return f
    
def glassify(snapshot, shape, changaPreset='default', verbose=False, \
             fulloutput=False):
    """
    Glassifies a snapshot, saves the results to the default filename (see 
    sphglass.filenames()) and returns the snapshot.  snapshot can be a filename
    or a pynbody SimSnap
    """
    
    inFile, fPrefix = filenames()
    paramname = fPrefix + '.param'
    
    if not isinstance(snapshot, str):
        
        snapshot.write(filename=inFile, fmt=pynbody.tipsy.TipsySnap)
        snapshotName = inFile
        
    else:
        
        snapshotName = snapshot
        snapshot = pynbody.load(snapshotName)
        
    try:
        
        param = makeParam(snapshot, shape, fulloutput)
        diskpy.utils.configsave(param, paramname, 'param')
        shutil.move(snapshotName, inFile)
        glass = reglassify(changaPreset, verbose, fulloutput)
        
    finally:
        
        shutil.move(inFile, snapshotName)
        
    return glass
    
def reglassify(changaPreset='default', verbose=True, fulloutput=False):
    """
    Run the most recently created glass (in the current working directory) 
    again to make it more glassy.
    
    Parameters
    ----------
    changaPreset : str
        ChaNGa preset to use (see diskpy for info on the ChaNGa presets)
    verbose : bool
        Verbosity, true or false
    fulloutput : bool
        If true, don't clean-up extra snapshot stuff
    
    Returns
    -------
    f : SimSnap (see pynbody)
        Snapshot containing the generated particle positions.  The glass
        positions can be accessed using f['pos'].  Density is also available for
        all the particles in f['rho'].  The target density is 1
    """
    
    # Get the default paramname
    paramname = filenames()[1] + '.param'
    # Glassify
    f = runchanga(paramname, changaPreset, verbose, fulloutput)
    return f
    
def runchanga(paramname, changaPreset='default', verbose=True, fulloutput=False):
    """
    Time evolves a snapshot in ChaNGa and overwrites the ICs with the result.
    Also sets the velocity to zero.
    
    Parameters
    ----------
    paramname : str
        Path to the .param file to run ChaNGa on
    changaPreset : str
        ChaNGa preset to use (see diskpy for info on the ChaNGa presets)
    verbose : bool
        Verbosity, true or false
    fulloutput : bool
        If true, don't clean-up extra snapshot stuff
    
    Returns
    -------
    f : SimSnap
        Simulation snapshot which has been run
    """
    
    param = diskpy.utils.configparser(paramname,'param')
    command = changa_command(paramname, changaPreset)
    p = changa_run(command, verbose=False, force_wait=False)
    
    currentStep = 0
    
    for line in iter(p.stdout.readline, ''):
            
            if verbose:
                print line,
            elif line[0:5] == 'Step:':
                
                line = line.strip()
                i = int(float(line.split()[1]))
                
                if i > currentStep:
                    
                    currentStep = i
                    print line, ' Total steps: ', param['nSteps']
    
    # move results and clean up
    fname = param['achOutName'] + '.{0:06}'.format(param['nSteps'])
    ICname = param['achInFile']
    if fulloutput:
        
        shutil.copyfile(fname, ICname)
        
    else:
        
        shutil.move(fname, ICname)
        os.system('rm -f ' + fname + '*')
        os.system('rm -f ' + param['achOutName'] + '.0*')
    
    # set velocity to zero
    f = pynbody.load(ICname, paramname=paramname)
    f['vel'] *= 0
    f.write()    
    
    return f
    
def boxSnap(n, shape, usegrid=False, randomness=0.):
    """
    Initialize snap shot with n randomly placed gas particles inside a box.
    
    if usegrid=True, a grid is used to seed the particle positions (see grid())
    n can then be [nx, ny, nz, ...] to specify the resolution along each
    dimension.  Otherwise, n is an int and a uniform spacing is attempted.
    """
    nDim = len(shape)
    if (nDim > 3) or (nDim < 1):
        
        raise ValueError, 'Only supported dimensions are 1, 2, 3. try'\
        'different shape'
        
    if usegrid:
        
        if hasattr(n, '__iter__'):
            
            res = n
            n = np.product(res)
            
        else:
            
            alpha = (float(n)/np.product(shape))**(1.0/nDim)
            res = np.array([alpha * L for L in shape])
            res = np.round(res).astype(int)
            n = np.product(res)
        
    n = int(n)
    snap = pynbody.new(gas=n)
    if usegrid:
        pos = grid(res, shape, randomness)
    else:
        pos = randomBox(n, shape)
    
    i0 = 3-nDim
    snap['pos'][:, i0:] = SimArray(pos,'au')
    volume = float(np.prod(shape))
    snap['mass'] = volume*SimArray(np.ones(n), 'Msol')/n
    snap['vel'] = SimArray(np.zeros([n,3]), 'km s**-1')
    snap['temp'] = SimArray(np.ones(n),'K')
    snap['eps'] = SimArray(np.ones(n))*max(shape)
    snap['rho'] = SimArray(np.ones(n), 'Msol kpc**-3')
    
    return snap
    
def getcs(snap, param):
    """
    From a simulation snapshot and param file (or dict), return the average
    sound speed. (in simulation units)
    """
    
    # Get sound speed
    units = diskpy.pychanga.units_from_param(param)
    mu = units['m_unit']
    tu = units['t_unit']
    lu = units['l_unit']
    m = SimArray(param['dMeanMolWeight'], 'm_p', dtype=float)
    m.convert_units(mu)
    K = pynbody.units.Unit('K')
    kBunit = mu*lu**2/(K*tu**2)
    k = kB.in_units(kBunit)
    T = snap['temp'].mean()
    cs2 = k*T/m
    cs = SimArray(float(np.sqrt(cs2)), lu/tu)
    
    return cs

def runTime(snap, param, boxShape):
    """
    From a SimSnap and a param file (or dict), estimate run time (in simulation
    units) required to evolve to a glass
    """
    h = estSmoothLength(snap, boxShape, param)
    cs = getcs(snap, param)
    t = 20 * h/cs
    return t

def estSmoothLength(snap, boxShape, param):
    """
    Get an estimate of a reasonable smoothing length for a snapshot
    """
    V = np.prod(boxShape)
    N = len(snap)
    nDim = len(boxShape)
    nSmooth = param.get('nSmooth', 32)
    h = 0.5 * ( 3*nSmooth*V/(4*np.pi*N))**(1./nDim)
    return h

def makeParam(snap, boxShape, fulloutput=False, runTimeScale=1., dampingScale=1.):
    """
    Make a param dict for creating a glass
    
    The number of dimensions, len(boxShape), must be 1, 2, or 3.
    """    
    # Get default params
    param = diskpy.utils.configparser(defaultparam)
    # Set Box size
    nDim = len(boxShape)
    maxL = max(boxShape)
    if nDim < 3:
        param['dxPeriod'] = 100 * maxL
    else:
        param['dxPeriod'] = float(boxShape[0])
    if nDim < 2:
        param['dyPeriod'] = 100 * maxL
    else:
        param['dyPeriod'] = float(boxShape[-2])
        
    param['dzPeriod'] = float(boxShape[-1])
    # Calculate run-time
    t = runTime(snap, param, boxShape) * runTimeScale
    param['dDelta'] = t/param['nSteps']
    # Setup output interval
    if fulloutput:
        
        param['iOutInterval'] = 1
        
    else:
        
        param['iOutInterval'] = param['nSteps']
        
    h = estSmoothLength(snap, boxShape, param)
    cs = getcs(snap, param)
    param['dGlassDamper'] = float(cs/h) * dampingScale
    
    return param
    
def randomBox(n, shape=[1,1,1]):
    """
    Generate random particle positions in an N-dimensional box of a given shape,
    cenetered at the origin
    """
    
    # Generate randomly placed particles inside the box
    nDim = len(shape)
    x = np.random.rand(n, nDim) - 0.5
    
    for i, L in enumerate(shape):
        
        x[:,i] *= L
    
    return x

def grid(res=[10, 10, 10], shape=[1,1,1], randomness=0., putInBox=True):
    """
    Creates a grid of positions shifted randomly in ND by a fixed amount 
    controlled by randomness.  The positions are shifted randomly (in each dim)
    by randomness * dx
    """
    if len(res) != len(shape):
        
        raise ValueError, 'res and shape length do not match'
        
    bounds = []
    dxs = []
    
    for L, n in zip(shape, res):
        
        x = np.linspace(-L/2., L/2., n+1)[1:]
        dx = x[1] - x[0]
        x -= dx/2.
        bounds.append(x)
        dxs.append(dx)
        
    pos = np.array([a for a in itertools.product(*bounds)])
    
    for i in range(len(res)):
        
        pos[:,i] += dxs[i] * randomness * (np.random.random(len(pos)) - 0.5)
    
    if putInBox:
        
        placeInBox(pos, shape)
        
    return pos
    
def placeInBox(pos, boxshape=[1,1,1]):
    """
    Anything outside of the box bounds specified by boxshape will be shifted
    by one box length.  This means anything more than L away from the box
    boundary will not be shifted to inside the box.
    
    IN PLACE
    """
    N, nDim = pos.shape       
    for i in range(nDim):
        
        L = boxshape[i]
        mask = pos[:, i] > L/2.
        pos[mask, i] -= L
        mask = pos[:, i] < -L/2.
        pos[mask, i] += L
        
