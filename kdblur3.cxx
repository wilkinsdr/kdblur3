// Subroutine to smooth the model spectrum by relativistic effects from a
// disk in the presence of a non-spinning black hole - uses laor3.
// ACF and RMJ May/June 1998
// kaa converted to C++ Aug 2018
// kaa changed to use FFTs Aug 2019
// DRW version for twice-broken power law emissivity, Aug 2019

//  parameters :
//       0        power law index for emissivity (10 for disk)
//       1        inner radius (GM/c**2)
//       2        outer radius (GM/c**2)
//       3        inclination  (degrees)
//       4        inner break radius (GM/c**2)
//       5        middle power-law dependence
//       6        outer break radius (GM/c**2)
//       7        outer power-law dependence

#include <xsTypes.h>
#include <functionMap.h>
#include <XSUtil/Numerics/Convolution.h>

extern "C" void laor3(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
	   const string& initString);

extern "C" void kdblur3(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{

  RealArray laor3Params(9);

  // set the laor2 energy to the midpoint of the energy range (this is arbitrary)
 
  laor3Params[0] = (energyArray[energyArray.size()-1]+energyArray[0])/2.0;
  for (size_t i=0; i<8; i++) laor3Params[i+1] = params[i];

  Numerics::ConvolutionInLnSpace<laor3>(energyArray, laor3Params, laor3Params[0], 
				       spectrumNumber, "", flux, fluxErr);

  return;
}
