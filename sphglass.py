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
import sys
import shutil

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
    
    return diskpy.utils.configparser(defaultparam, 'param')
    
def filenames():
    
    param = _loadDefaults()
    inFile = param['achInFile']
    outPrefix = param['achOutName']
    
    return inFile, outPrefix

def glassBox(n, shape=[1,1,1], changaPreset='default', verbose=False, 
              fulloutput=False, accuracy=0.25, max_reglass=50):
    """
    Generates an sph glass in a box with periodic boundary conditions using 
    ChaNGa.  The procedure is:
        1. Generate random particle positions in a box
        2. Create a tipsy snapshot with only gas particles
        3. Time evolve in a periodic box with no gravity and lots of artificial
            viscosity.
    
    Parameters
    ----------
    n : int
        Number of particles
    shape : array-like
        Shape of the box (x, y, z).  The box will be centered around the origin.
    changaPreset : str
        ChaNGa preset to use (see diskpy for info on the ChaNGa presets)
    verbose : bool
        Verbosity, true or false
    fulloutput : bool
        If True, all the snapshots for each time step during the time evolution
        will be output.
    accuracy : float
        A parameter to scale the required accuracy for glass convergence.  The
        glass is evoled until (stddev(rho) * nSmooth) < accuracy
    max_reglass : int
        The maximum number of times to re-run the time evolution using reglassify().
        Each time the snapshot is time evolved the velocities are set to zero.
        This is useful especially with long boxes where waves can form.
    
    Returns
    -------
    f : SimSnap (see pynbody)
        Snapshot containing the generated particle positions.  The glass
        positions can be accessed using f['pos'].  Density is also available for
        all the particles in f['rho'].  The target density is 1
    
    Notes
    -----
    If you find the box is not sufficiently glassy, you can time evolve it
    again by running reglassify()
    """
    if max_reglass is None:
        
        max_reglass = np.inf
        
    # Generate snapshot with random positions
    snap = boxSnap(n, shape)
    param = makeParam(snap, shape, fulloutput)
    # Save snapshot and param
    paramname = param['achOutName'] + '.param'
    ICname = param['achInFile']
    nSmooth = param['nSmooth']
    diskpy.utils.configsave(param, paramname)
    snap.write(fmt=pynbody.tipsy.TipsySnap, filename=ICname)
    # Run ChaNGa to make a glass
    f = runchanga(paramname, changaPreset, verbose, fulloutput)
    
    i = 1
    
    while i < max_reglass:
        
        relvar = f['rho'].std() * nSmooth
        print 'ITERATION:', i
        print 'TARGET ACCURACY:', accuracy
        print 'CURRENT ACCURACY:', relvar
        sys.stdout.flush()
        
        if relvar <= accuracy:
            
            # Success, we have converged to a glass
            print 'Successfully converged to a glass in {0} iterations'.format(i)
            break
            
        f = reglassify(changaPreset, verbose, fulloutput)
        i += 1
        
    
#    for i in range(nreglass):
#        
#        f = reglassify(changaPreset, verbose, fulloutput)
    
    return f
    
def glassify(snapshot, shape, changaPreset='default', verbose=False, fulloutput=False):
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
    
def boxSnap(n, shape):
    """
    Initialize snap shot with n randomly placed gas particles inside a box.
    """
    n = int(n)
    snap = pynbody.new(gas=n)
    snap['pos'] = SimArray(randomBox(n, shape),'au')
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

def runTime(snap, param):
    """
    From a SimSnap and a param file (or dict), estimate run time (in simulation
    units) required to evolve to a glass
    """
    cs = getcs(snap, param)
    shape = snap['pos'].max(0) - snap['pos'].min(0)
    particleVol = np.prod(shape)/len(snap)
    L = particleVol**(1,3)
    
    t = float(64 * L/cs)
    # multiply by a 'fudge' factor
    t *= 0.3
    return t

def makeParam(snap, boxShape, fulloutput=False):
    """
    Make a param dict for creating a glass
    """    
    # Get default params
    param = diskpy.utils.configparser(defaultparam)
    # Set Box size
    param['dxPeriod'] = float(boxShape[0])
    param['dyPeriod'] = float(boxShape[1])
    param['dzPeriod'] = float(boxShape[2])
    # Calculate run-time
    t = runTime(snap, param)
    param['dDelta'] = t/param['nSteps']
    
    # Setup output interval
    if fulloutput:
        
        param['iOutInterval'] = 1
        
    else:
        
        param['iOutInterval'] = param['nSteps']
    
    # Check that the courant condition isn't too weak.  A very weak courant
    # condition can result in particles traveling more than 1 box-length
    # in a single step.  This raises an error in ChaNGa
    V = np.prod(boxShape)
    N = len(snap)
    nSmooth = param.get('nSmooth', 32)
    h = 0.5 * ( 3*nSmooth*V/(4*np.pi*N))**(1./3)
    h = float(h) * (1 + 1./np.sqrt(N)) # include correction for poisson noise
    L = min(boxShape)
    
    courantMax = L/(50*h)
    
    if param['dEtaCourant'] > courantMax:
        
        param['dEtaCourant'] = courantMax
        print 'Low resolution.  Dropping dEtaCourant to ', courantMax
        print 'This is to avoid particles moving more than one box length in a tstep'
        
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
    
def _rhovar():
    """
    gets stddev(rho)/mean(rho) vs time
    
    returns rhostd, t
    """
    
    fnames = diskpy.pychanga.get_fnames('glass')
    fs = []
    for fname in fnames:
        fs.append(pynbody.load(fname, paramname='glass.param'))
    rhos = np.zeros(len(fs))    
    t = np.zeros(len(fs))
    for i, f in enumerate(fs):
        
        rhos[i] = f['rho'].std()/f['rho'].mean()
        t[i] = diskpy.pychanga.snapshot_time(f)
    
    return rhos, t