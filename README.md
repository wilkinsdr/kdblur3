# kdblur3
XSPEC models for modelling the relativistically broadened emission lines and blurred reflection arising from an accretion disc illuminated according to a twice-broken power law emissivity profile, described in Wilkins & Fabian 2011, MNRAS 414, 1269.

laor3 is the emission line model. kdblur3 is the convolution kernel for blurring a rest-frame reflection spectrum.

The XSPEC model function is written in C++ (compconv.cxx) but calls Fortran routines to perform the calculation (docompconf.f). The XSPEC model definition is in lmodel_kdblur3.dat.

Parameters
----------
1. lineE   : (laor3 only) Rest frame energy of the emission line (keV)
2. Index   : Power law index of the inner section of the emissivity profile (emissivity goes as r^-q)
3. Rin(G)  : Inner radius of the accretion disc (gravitaional radii)
4. Rout(G) : Outer radius of the accretion disc (gravitaional radii)
5. Incl    : Inclination of the observer's line of sight to the disc normal (deg)
6. Rbreak1 : Inner break radius of twice-broken power law emissivity profile
7. Index1  : Power law index of the middle section of the emissivity profile
8. Rbreak2 : Outer break radius of twice-broken power law emissivity profile
9. Index2  : Power law index of the outer section of the emissivity profile

Installation
------------
Installation

The model is compiled using the initpackage command in XSPEC:

XSPEC> cd kdblur3/model/directory

XSPEC> initpackage kdblurthree lmodel_kdblur3.dat .

Then can be loaded from the directory it is compiled in (it will need to be loaded every time XSPEC is started):

XSPEC> lmod kdblurthree .

Alternatively, the code and model definition can be integrated into a local model package.
