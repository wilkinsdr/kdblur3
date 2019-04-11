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
void blurring(const RealArray& kernelEnergyArray, const RealArray& kernelFluxArray,
	      const Real kernelLineEnergy, const RealArray& energyArray,
	      RealArray& fluxArray, const bool Reset);

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

  // set the laor3 energy to 10keV (this is arbitrary)
 
  laor3Params[0] = 10.0;
  for (size_t i=0; i<8; i++) laor3Params[i+1] = params[i];

  // find the minimum binsize in the input energyArray
  Real minBinsize = energyArray[energyArray.size()-1] - energyArray[0];
  for (size_t i=0; i<energyArray.size()-1; i++) {
    Real binsize = energyArray[i+1] - energyArray[i];
    if ( binsize > 0.0 && binsize < minBinsize ) minBinsize = binsize;
  }

  // set up the energy array on which the diskline is calculated
  // from 0 to 20 keV in linear bins of size minBinsize.

  size_t kbins = (size_t)(20.0/minBinsize);
  RealArray kernelEnergyArray(kbins+1);
  for (size_t i=0; i<=kbins; i++) kernelEnergyArray[i] = (20.0*i)/((Real)kbins);

  laor3(kernelEnergyArray, laor3Params, spectrumNumber, kernelFluxArray,
       kernelFluxErrArray, initString);

  blurring(kernelEnergyArray, kernelFluxArray, laor3Params[0], energyArray, flux, true);

  return;
  
}
