/*
 
 laor3.cxx
 
 XSPEC model for relativistically broadened emission line from an accretion 
 disc with a twice-broken power law emissivity profile.
 
 Based on laor2 by A.C. Fabian and R.M. Johnstone.
 
 v2.0, November 2014, D.R. Wilkins
 
 */
#include <xsTypes.h>
#include <stlToCArrays.h>
#include <functionMap.h>
#include <XSUtil/FunctionUtils/FunctionUtility.h>
#include <XSUtil/FunctionUtils/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <fitsio.h>
#include <sstream>
#include <cfortran.h>

PROTOCCALLSFSUB17(LAORSB3,laorsb3,FLOATV,INT,FLOATV,INT,INT,INT,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV)

#define LAORSB3(ear,ne,param,nenrgy,nrad,nteta,efil,ebin,rad,incl,itrs,bin,start,end,fstart,fend,photar) \
           CCALLSFSUB17(LAORSB3,laorsb3,FLOATV,INT,FLOATV,INT,INT,INT,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV, \
           ear,ne,param,nenrgy,nrad,nteta,efil,ebin,rad,incl,itrs,bin,start,end,fstart,fend,photar)

namespace {
   void errReturn(fitsfile* pFile, const string& context, int status);
}


extern "C" void laor3(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{

   // This subroutine calculates line profiles from an accretion disk around
   // a Kerr BH with a twice-broken power law emissivity profile.

   // Parameters :
   //    1:    Rest frame energy of line
   //    2:    Power-law dependence of emission
   //    3:    Inner radius of disk (in gravl radii)  >= 1.235
   //    4:    Outer radius of disk (in gravl radii)  =< 400
   //    5:    Inclination angle of disk (degrees, 0=face on)
   //    6:    Break radius
   //    7:    outer power-law dependence

   // The transfer function is read from the file ari.bin
   // The Radii are in the range 1.23 to 400, inclinations are from 0 to 1.
   // There are 31 tabulated inclinations in steps of 0.033 from 0 to 1.

   // kaa 1/17/95 modified to use FITS file.
   //  kaa 9/6/96  added dynamic memory

   // Tabulated data:
   //    nrad         Number of radii
   //    nteta        Number of inclinations
   //    NENRGY       Number of energies
   //    radpt        Radii (in gravl radii)
   //    inclpt       cos(inclination)
   //    ebin         Energies (units of rest energy of line)
   //    trspt        Transfer function

   // Translated from laor2.f Jan 2009 (CG)

   using namespace XSutility;

   static auto_array_ptr<float> apRad(0);
   static auto_array_ptr<float> apIncl(0);
   static auto_array_ptr<float> apEbin(0);
   static auto_array_ptr<float> apEfil(0);
   static auto_array_ptr<float> apTrs(0);
   static auto_array_ptr<float> apBin(0);

   static bool isFirst = true;
   static int nrad=0;
   static int nenergy=0;
   static int nteta=0;

   // If the first time round then read in the data 

   if (isFirst)
   {
      string laorFile(FunctionUtility::modelDataPath());
      const string fName("ari.mod");
      laorFile += "ari.mod";
      fitsfile *pFile=0;
      int status=0;
      if (fits_open_file(&pFile, const_cast<char*>(laorFile.c_str()), 0, &status))
      {
         std::ostringstream err;
         err << "Failed to open " << laorFile << "\nFITSIO error = "
            << status << " in LAOR3\n";
         xs_write(const_cast<char*>(err.str().c_str()), 10);
         return;
      }

      // Go to the second extension to get the radii

      if (fits_movabs_hdu(pFile, 2, 0, &status))
      {
         string context("Failed to open second extension in " + fName);
         errReturn(pFile, context, status);
         return;
      }

      long lnrad=0;
      if (fits_read_key_lng(pFile, "NAXIS2", &lnrad, 0, &status))
      {
         string context("Failed to get size of RADIUS data");
         errReturn(pFile, context, status);
         return;
      }

      int anynul=0;
      float* oldF = apRad.reset(new float[lnrad]);
      delete [] oldF;
      if (fits_read_col_flt(pFile, 1, 1, 1, lnrad, 0., apRad.get(),
                &anynul, &status))
      {
         string context("Failed to read RADIUS data");
         errReturn(pFile, context, status);
         return;
      }

      // Go to the third extension to get the angles

      if (fits_movabs_hdu(pFile, 3, 0, &status))
      {
         string context("Failed to open third extension in " + fName);
         errReturn(pFile, context, status);
         return;
      }

      long lnteta=0; 
      if (fits_read_key_lng(pFile, "NAXIS2", &lnteta, 0, &status))
      {
         string context("Failed to get size of ANGLE data");
         errReturn(pFile, context, status);
         return;
      }

      oldF = apIncl.reset(new float[lnteta]);
      delete [] oldF;
      if (fits_read_col_flt(pFile, 1, 1, 1, lnteta, 0., apIncl.get(),
                &anynul, &status))
      {
         string context("Failed to read ANGLE data");
         errReturn(pFile, context, status);
         return;
      }

      // Go to the fourth extension to get the energies

      if (fits_movabs_hdu(pFile, 4, 0, &status))
      {
         string context("Failed to open fourth extension in " + fName);
         errReturn(pFile, context, status);
         return;
      }

      long lnenergy=0;
      if (fits_read_key_lng(pFile, "NAXIS2", &lnenergy, 0, &status))
      {
         string context("Failed to get size of ENERGY data");
         errReturn(pFile, context, status);
         return;
      }

      oldF = apEbin.reset(new float[lnenergy]);
      delete [] oldF;
      if (fits_read_col_flt(pFile, 1, 1, 1, lnenergy, 0., apEbin.get(),
                &anynul, &status))
      {
         string context("Failed to read ENERGY data");
         errReturn(pFile, context, status);
         return;
      }

      // Get the memory for the efilpt array
      oldF = apEfil.reset(new float[lnenergy+3]);
      delete [] oldF;

      // Go to the fifth extension to read the transfer function

      if (fits_movabs_hdu(pFile, 5, 0, &status))
      {
         string context("Failed to open fifth extension in " + fName);
         errReturn(pFile, context, status);
         return;
      }

      oldF = apTrs.reset(new float[lnrad*lnteta*lnenergy]);
      delete [] oldF;
      if (fits_read_col_flt(pFile, 1, 1, 1, lnrad*lnteta*lnenergy, 0., 
                apTrs.get(), &anynul, &status))
      {
         string context("Failed to read transfer function");
         errReturn(pFile, context, status);
         return;
      }

      fits_close_file(pFile, &status);

      oldF = apBin.reset(new float[lnenergy+2]);
      delete [] oldF;

      nrad = static_cast<int>(lnrad);
      nteta = static_cast<int>(lnteta);
      nenergy = static_cast<int>(lnenergy);

   } // end if isFirst

   // allocate memory for inibin/erebin work arrays
   int nE = static_cast<int>(energyArray.size()) - 1;
   auto_array_ptr<float> apStart(new float[nE]);
   auto_array_ptr<float> apEnd(new float[nE]);
   auto_array_ptr<float> apFstart(new float[nE]);
   auto_array_ptr<float> apFend(new float[nE]);

   // call the subroutine to do the actual calculation

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   auto_array_ptr<float> apEar(ear);
   auto_array_ptr<float> apPars(pars);
   auto_array_ptr<float> apPhotar(photar);
   auto_array_ptr<float> apPhoter(photer);

   LAORSB3(ear, nE, pars, nenergy, nrad, nteta, apEfil.get(), apEbin.get(),
        apRad.get(), apIncl.get(), apTrs.get(), apBin.get(), apStart.get(),
        apEnd.get(), apFstart.get(), apFend.get(), photar);

   XSFunctions::floatFluxToStl<float>(photar, photer, nE, false, flux, fluxErr);       

   isFirst = false;      
}

namespace {
   void errReturn(fitsfile* pFile, const string& context, int status)
   {
      std::ostringstream oss;
      oss << context << "\nError in LAOR3 : status = " << status;
      string fullMsg(oss.str());
      xs_write(const_cast<char*>(fullMsg.c_str()), 10);
      int status2=0;
      fits_close_file(pFile, &status2);
   }
}



void cppModelWrapper(const double* energy, int nFlux, const double* params,
					 int spectrumNumber, double* flux, double* fluxError, const char* initStr,
					 int nPar, void (*cppFunc)(const RealArray&, const RealArray&,
											   int, RealArray&, RealArray&, const string&));

void fcppModelWrapper(const float* energy, int nFlux, const float* params,
					  int spectrumNumber, float* flux, float* fluxError,
					  int nPar, void (*cppFunc)(const RealArray&, const RealArray&,
												int, RealArray&, RealArray&, const string&));


void cppModelWrapper(const double* energy, int nFlux, const double* params,
					 int spectrumNumber, double* flux, double* fluxError, const char* initStr,
					 int nPar, void (*cppFunc)(const RealArray&, const RealArray&,
											   int, RealArray&, RealArray&, const string&))
{
	// Assumes energy points to arrays of size nFlux+1, flux and fluxError
	// point to arrays of size nFlux (though they need not be initialized),
	// and params points to an array of size nPar.
	RealArray energy_C(energy, (size_t)nFlux+1);
	RealArray params_C(params, nPar);
	RealArray flux_C(flux, (size_t)nFlux);
	RealArray fluxError_C(fluxError, (size_t)nFlux);
	string cppStr;
	if(initStr && strlen(initStr))
	cppStr = initStr;
	(*cppFunc)(energy_C, params_C, spectrumNumber, flux_C, fluxError_C, cppStr);
	for (int i=0; i<nFlux; ++i)
	{
		flux[i] = flux_C[i];
	}
	if (fluxError_C.size())
	{
		for (int i=0; i<nFlux; ++i)
		{
			fluxError[i] = fluxError_C[i];
		}
	}
}

void fcppModelWrapper(const float* energy, int nFlux, const float* params,
					  int spectrumNumber, float* flux, float* fluxError,
					  int nPar, void (*cppFunc)(const RealArray&, const RealArray&,
												int, RealArray&, RealArray&, const string&))
{
	// Assumes energy points to arrays of size nFlux+1, flux and fluxError
	// point to arrays of size nFlux (though they need not be initialized),
	// and params points to an array of size nPar.
	RealArray energy_C(0.0, (size_t)nFlux+1);
	RealArray params_C(0.0, nPar);
	RealArray flux_C(0.0, (size_t)nFlux);
	RealArray fluxError_C(0.0, (size_t)nFlux);
	string cppStr;
	
	for (int i=0; i<nFlux+1; ++i)
	{
		energy_C[i] = (double)energy[i];
	}
	for (int i=0; i<nPar; ++i)
	{
		params_C[i] = (double)params[i];
	}
	for (int i=0; i<nFlux; ++i)
	{
		flux_C[i] = (double)flux[i];
	}
	if (fluxError)
	{
		for (int i=0; i<nFlux; ++i)
		{
			fluxError_C[i] = (double)fluxError[i];
		}
	}
	
	
	(*cppFunc)(energy_C, params_C, spectrumNumber, flux_C, fluxError_C, cppStr);
	
	for (int i=0; i<nFlux; ++i)
	{
		flux[i] = (float)flux_C[i];
	}
	if (fluxError_C.size())
	{
		for (int i=0; i<nFlux; ++i)
		{
			fluxError[i] = (float)fluxError_C[i];
		}
	}
}

void f_laor3(const float* energy, int nFlux, const float* params,
			 int spectrumNumber, float* flux, float* fluxError);

void f_laor3(const float* energy, int nFlux, const float* params,
			 int spectrumNumber, float* flux, float* fluxError)
{
	const size_t nPar = 9;
	fcppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
					 nPar, laor3);
}

void C_laor3(const double* energy, int nFlux, const double* params,
			 int spectrumNumber, double* flux, double* fluxError, const char* initStr)
{
	const size_t nPar = 9;
	cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
					initStr, nPar, laor3);
}

FCALLSCSUB6(f_laor3,LAOR3,laor3,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
FCALLSCSUB7(C_laor3,DLAOR3,dlaor3,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,DOUBLEV,STRING)



