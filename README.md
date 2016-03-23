# sphglass
A simple package for generating an SPH glass in a periodic box written by Isaac Backus.

# Installing sphglass
Download and install depedencies:
 * [ChaNGa](https://github.com/N-BodyShop/changa)
 * [pynbody](https://github.com/pynbody/pynbody)
 * [diskpy](https://github.com/ibackus/diskpy)

Download:

`git clone https://github.com/ibackus/sphglass.git`

Then add the `sphglass/` directory to your `PYTHONPATH`.

# Example

```python
import sphglass
import matplotlib.pyplot as plt

# Generate a glass in a box shape 1 by 2 by 3 (in x,y,z).  f is a SimSnap
# (see pynbody).  The only relevant quantites are the position and density
f = sphglass.glassBox(1000, [1,2,3])
# A good glass should have nearly uniform density, although there
# will always be some scatter
plt.plot(f['z'], f['rho'], 'x')
# If you think this is not uniform enough, you try to make it more glassy:
f = sphglass.reglassify()
plt.plot(f['z'], f['rho'], 'x')
# Otherwise, you can just re-glassify many times (re-run the time evolution)
f = sphglass.glassbox(1000, [1,2,3], nreglass=5)
```

# More information
This works best on boxes that are close to cubic.  Waves can be generated in long, narrow boces.  The algorithm used here for generating a glass is pretty simple. The procedure is:
 1. Generate random particle positions in a box
 2. Create a tipsy snapshot with only gas particles
 3. Time evolve in a periodic box with no gravity and lots of artificial
  viscosity.
A random distribution (Poisson distribution) of particles is in a higher energy state (higher density) and will gradually relax to a lower energy, glass-like, uniform density state given enough time.  The trick is to just estimate how much time is required for that and to have a dissipative force, in this case lots of artificial viscosity.

The default parameters used for this are stored in glassdefaults.param and can be edited.
