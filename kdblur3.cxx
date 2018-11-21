// Subroutine to smooth the model spectrum by relativistic effects from a
// disk in the presence of a non-spinning black hole - uses laor3.
// ACF and RMJ May/June 1998
// kaa converted to C++ Aug 2018
// DRW modified kdblur2->kdblur3 October 21, 2018

//  parameters :
//       0        power law index for emissivity (10 for disk)
//       1        inner radius (GM/c**2)
//       2        outer radius (GM/c**2)
//       3        inclination  (degrees)
//       4        break radius (GM/c**2)
//       5        outer power-law dependence

#include <xsTypes.h>
#include <functionMap.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSstreams.h>

// from blurring.cxx
void blurring(const RealArray& energyArray, const RealArray& fluxArrayKernel,
	      const Real kernelLineEnergy, RealArray& fluxArray, const bool Reset);

// from laor3.cxx
extern "C" void laor3(const RealArray& energyArray, const RealArray& params,
           int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
           const string& initString);


extern "C" void kdblur3(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{

  RealArray kernelFluxArray, kernelFluxErrArray;

  RealArray laor3Params(9);

  // set the laor2 energy to about the middle of the input energy range
 
  laor3Params[0] = energyArray[energyArray.size()/2];
  for (size_t i=0; i<8; i++) laor3Params[i+1] = params[i];

  laor3(energyArray, laor3Params, spectrumNumber, kernelFluxArray,
       kernelFluxErrArray, initString);

  blurring(energyArray, kernelFluxArray, laor3Params[0], flux, true);

  return;
  
}
