/**
 * @file acquisition.cpp
 *
 * @brief acquisition source file for SDR functions in C++
 *
 * Project Title: GNSS-R SDR
 * Author: John Bagshaw
 * Co-author: Surabhi Guruprasad
 * Contact: jotshaw@yorku.ca
 * Supervisor: Prof. Sunil Bisnath
 * Project Manager: Junchan Lee
 * Institution: Lassonde School of Engineering, York University, Canada.
 **/


 /***********************************
 * Includes
 ***********************************/
#include "acquisition.hpp"
#include <math.h>
#include "gen_code.hpp"
#include <string.h>
#include "constants.h"
#include "fftw3.h"
#include "result_read_write.hpp"
#include "timing.hpp"
#include <thread>
#include <mutex>

 /***********************************
 * Defines
 ***********************************/

 /***********************************
 * Macros
 ***********************************/
std::mutex mtx;

void fft(
  fftw_plan p
) 
{ 
  fftw_execute(p);
}


 /***********************************
 * Static Variables
 ***********************************/

 /***********************************
 * Static Function Definitions
 ***********************************/
static void subtractMean(double* y, const double* x, const double len)
{
  double accum;
  double mean;
  int k;
  accum = 0.0;
  for (k = 0; k < len; k++) {
    accum += x[k];
  }
  mean = accum / len;

  for (k = 0; k < len; k++) {
    y[k] = x[k] - mean;
  }
}

void acquisitionCore(
  const Settings_t* settings,
  AcqResults_t* acqResults,
  int    PRNIdx,
  double* signal0DC,
  const double* longSignal,
  double** caCodeTables,
  bool* sigFound
)
{
#define FUNC_PREFIX " Normal Acquisition (): "
  mtx.lock();
  disp("PRN: %i\n", PRNIdx + 1);
  mtx.unlock();

  // -- - Initialize arrays to speed up the code------------------------------ -
  // Search results of all frequency bins and code shifts(for one satellite)

  // Number of the frequency bins for the given acquisition band(500Hz steps)
  double samplesPerCode = round(settings->samplingFreq / (settings->codeFreqBasis / settings->codeLength));
  int numberOfFrqBins   = round(settings->acqSearchBand * 2) + 1;
  double tS             = 1 / settings->samplingFreq;
  double tC             = 1.0 / settings->codeFreqBasis;
  int fftNumPts         = (int)(8 * pow(2, ceil(log(10 * samplesPerCode) / log(2))));
  int uniqFftPts        = (int)ceil((fftNumPts + 1) / 2.0);

  double* xCarrier    = new double  [10 * samplesPerCode];
  double* fftFreqBins = new double  [uniqFftPts];
  double** results    = new double* [numberOfFrqBins];

  double acqMS = floor(settings->acqMS / 2.0);

  fftw_complex* caCodeFreqDomIn = new fftw_complex   [samplesPerCode];
  fftw_complex* caCodeFreqDom   = new fftw_complex   [samplesPerCode];
  fftw_complex* acqRes1Temp     = new fftw_complex   [samplesPerCode];
  fftw_complex* acqRes2Temp     = new fftw_complex   [samplesPerCode];
  fftw_complex* acqRes1         = new fftw_complex   [samplesPerCode];
  fftw_complex* acqRes2         = new fftw_complex   [samplesPerCode];
  fftw_complex* xCarrierIn      = new fftw_complex   [fftNumPts     ];
  fftw_complex* xCarrierOut     = new fftw_complex   [fftNumPts     ];
  fftw_complex*** IQfreqDom1In  = new fftw_complex **[acqMS];
  fftw_complex*** IQfreqDom2In  = new fftw_complex **[acqMS];
  fftw_complex*** IQfreqDom1    = new fftw_complex **[acqMS];
  fftw_complex*** IQfreqDom2    = new fftw_complex **[acqMS];
  fftw_complex* IQfreqDom1Loc   = new fftw_complex   [samplesPerCode];
  fftw_complex* IQfreqDom2Loc   = new fftw_complex   [samplesPerCode];

  fftw_plan p1 = fftw_plan_dft_1d((int)samplesPerCode, caCodeFreqDomIn, caCodeFreqDom, false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d((int)samplesPerCode, IQfreqDom1Loc  , acqRes1      , false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p3 = fftw_plan_dft_1d((int)samplesPerCode, IQfreqDom2Loc  , acqRes2      , false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p4 = fftw_plan_dft_1d(fftNumPts     , xCarrierIn     , xCarrierOut  , true  ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

  double scale = 1.0 / samplesPerCode;
  for (size_t m = 0; m < acqMS; m++)
  {
    double* signal1 = new double[samplesPerCode];
    double* signal2 = new double[samplesPerCode];

    memcpy(signal1, &longSignal[(2 * m + 0) * (size_t)samplesPerCode], sizeof(double) * samplesPerCode);
    memcpy(signal2, &longSignal[(2 * m + 1) * (size_t)samplesPerCode], sizeof(double) * samplesPerCode);

    IQfreqDom1In[m] = new fftw_complex * [numberOfFrqBins];
    IQfreqDom2In[m] = new fftw_complex * [numberOfFrqBins];
    IQfreqDom1  [m] = new fftw_complex * [numberOfFrqBins];
    IQfreqDom2  [m] = new fftw_complex * [numberOfFrqBins];

    for (int i = 0; i < numberOfFrqBins; ++i)
    {
      results[i] = new double[samplesPerCode];
      IQfreqDom1In[m][i] = new fftw_complex[samplesPerCode];
      IQfreqDom2In[m][i] = new fftw_complex[samplesPerCode];
      IQfreqDom1  [m][i] = new fftw_complex[samplesPerCode];
      IQfreqDom2  [m][i] = new fftw_complex[samplesPerCode];

      fftw_plan iqFreqDom1Plans = fftw_plan_dft_1d((int)samplesPerCode, &IQfreqDom1In[m][i][0], &IQfreqDom1[m][i][0], true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
      fftw_plan iqFreqDom2Plans = fftw_plan_dft_1d((int)samplesPerCode, &IQfreqDom2In[m][i][0], &IQfreqDom2[m][i][0], true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

      double frqBins = settings->iF - (settings->acqSearchBand / 2) * 1000 + 0.5e3 * (i);
      for (size_t j = 0; j < samplesPerCode; j++)
      {
        double sinCar = sin(frqBins * j * 2 * Pi * tS);
        double cosCar = cos(frqBins * j * 2 * Pi * tS);

        // Generate carrier wave frequency grid(0.5kHz step)-----------
        // -- - Generate local sine and cosine------------------------------
        // -- - Convert the baseband signal to frequency domain-------------

        IQfreqDom1In[m][i][j][0] = sinCar * signal1[j];
        IQfreqDom1In[m][i][j][1] = cosCar * signal1[j];
        IQfreqDom2In[m][i][j][0] = sinCar * signal2[j];
        IQfreqDom2In[m][i][j][1] = cosCar * signal2[j];

      }

      fft(iqFreqDom1Plans);
      fft(iqFreqDom2Plans);

      fftw_destroy_plan(iqFreqDom1Plans);
      fftw_destroy_plan(iqFreqDom2Plans);

    } // frequency bins

    delete[] signal1;
    delete[] signal2;

  }

  // -- - Initialize acqResults------------------------------------------------
  
  // Carrier frequencies of detected signals
  acqResults->carrFreq[PRNIdx] = 0.0;
  // C / A code phases of detected signals
  acqResults->codePhase[PRNIdx] = 0.0;
  // Correlation peak ratios of the detected signals
  acqResults->peakMetric[PRNIdx] = 0.0;

  int PRN = (int)settings->acqSatelliteList[PRNIdx];

    //// Correlate signals ======================================================
    for (size_t i = 0; i < (size_t)samplesPerCode; i++)
    {
      caCodeFreqDomIn[i][0] = caCodeTables[PRNIdx][i];
      caCodeFreqDomIn[i][1] = 0.0;
    }
    // TAKE IFFT HERE
    fft(p1);

    // Make the correlation for whole frequency band(for all freq.bins)
    for (size_t frqBinIndex = 0; frqBinIndex < numberOfFrqBins; frqBinIndex++)
    {
      //mtx.lock();
      //disp("PRN: %i / %i, frqBinIndex: %i / %i\n", PRNIdx + 1, 32, frqBinIndex + 1, numberOfFrqBins);
      //mtx.unlock();

      // -- - Multiplication in the frequency domain(correlation in time
      // domain)

      for (size_t j = 0; j < samplesPerCode; j++)
      {
        acqRes1Temp[j][0] = 0.0;
        acqRes1Temp[j][1] = 0.0;
        acqRes2Temp[j][0] = 0.0;
        acqRes2Temp[j][1] = 0.0;
      }

      for (size_t m = 0; m < acqMS; m++)
      {

        for (size_t i = 0; i < samplesPerCode; i++)
        {
          IQfreqDom1Loc[i][0] = (IQfreqDom1[m][frqBinIndex][i][0] * caCodeFreqDom[i][0]) - (IQfreqDom1[m][frqBinIndex][i][1] * caCodeFreqDom[i][1]);
          IQfreqDom1Loc[i][1] = (IQfreqDom1[m][frqBinIndex][i][0] * caCodeFreqDom[i][1]) + (IQfreqDom1[m][frqBinIndex][i][1] * caCodeFreqDom[i][0]);

          IQfreqDom2Loc[i][0] = (IQfreqDom2[m][frqBinIndex][i][0] * caCodeFreqDom[i][0]) - (IQfreqDom2[m][frqBinIndex][i][1] * caCodeFreqDom[i][1]);
          IQfreqDom2Loc[i][1] = (IQfreqDom2[m][frqBinIndex][i][0] * caCodeFreqDom[i][1]) + (IQfreqDom2[m][frqBinIndex][i][1] * caCodeFreqDom[i][0]);
        }
        // TAKE IFFT HERE
        fft(p2);
        fft(p3);
        for (size_t i = 0; i < samplesPerCode; i++)
        {
          acqRes1[i][0] *= scale;
          acqRes2[i][0] *= scale;
          acqRes1[i][1] *= scale;
          acqRes2[i][1] *= scale;

          acqRes1Temp[i][0] += acqRes1[i][0];
          acqRes2Temp[i][0] += acqRes2[i][0];
          acqRes1Temp[i][1] += acqRes1[i][1];
          acqRes2Temp[i][1] += acqRes2[i][1];
        }
      }

      // -- - Perform inverse DFT and store correlation results------------
      double maxAcqRes1 = 0.0;
      double maxAcqRes2 = 0.0;

      for (size_t i = 0; i < samplesPerCode; i++)
      {
        acqRes1Temp[i][0] = acqRes1Temp[i][0] * acqRes1Temp[i][0] + acqRes1Temp[i][1] * acqRes1Temp[i][1];
        acqRes2Temp[i][0] = acqRes2Temp[i][0] * acqRes2Temp[i][0] + acqRes2Temp[i][1] * acqRes2Temp[i][1];

        if (acqRes1Temp[i][0] > maxAcqRes1)
        {
          maxAcqRes1 = acqRes1Temp[i][0];
        }

        if (acqRes2Temp[i][0] > maxAcqRes2)
        {
          maxAcqRes2 = acqRes2Temp[i][0];
        }
      }

      for (size_t i = 0; i < samplesPerCode; i++)
      {
        if (maxAcqRes1 > maxAcqRes2)
        {
          results[frqBinIndex][i] = acqRes1Temp[i][0];
        }
        else
        {
          results[frqBinIndex][i] = acqRes2Temp[i][0];
        }
      }
    } // frequency bin


    // //Look for correlation peaks in the results ==============================
    // Find the highest peak and compare it to the second highest peak
    // The second peak is chosen not closer than 1 chip to the highest peak

    // -- - Find the correlation peak and the carrier frequency--------------
    // -- - Find code phase of the same correlation peak-------------------- -
    double maxPeakSize = 0.0;
    double peakSize = 0.0;
    int frequencyBinIndex = 0;
    int codePhase = 0;
    double* maxF = new double[numberOfFrqBins];

    for (size_t f = 0; f < numberOfFrqBins; f++)
    {
      maxPeakSize = 0.0;
      for (size_t i = 0; i < samplesPerCode; i++)
      {
        if (results[f][i] > maxPeakSize)
        {
          maxPeakSize = results[f][i];
        }
      }
      maxF[f] = maxPeakSize;
    }

    peakSize = 0.0;
    for (size_t f = 0; f < numberOfFrqBins; f++)
    {
      if (maxF[f] > peakSize)
      {
        peakSize = maxF[f];
        frequencyBinIndex = f;
      }
    }

    maxPeakSize = 0.0;
    for (size_t f = 0; f < numberOfFrqBins; f++)
    {
      for (size_t i = 0; i < samplesPerCode; i++)
      {
        if (results[f][i] > maxPeakSize)
        {
          maxPeakSize = results[f][i];
          codePhase = i;
        }
      }
    }

    delete[] maxF;

    // -- - Find 1 chip wide C / A code phase exclude range around the peak----
    double samplesPerCodeChip = round(settings->samplingFreq / settings->codeFreqBasis);
    double excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    double excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    // -- - Correct C / A code phase exclude range if the range includes array
    // boundaries
    int* codePhaseRange = new int[samplesPerCode];
    int newerCodeLen = 0;
    memset(codePhaseRange, 0, sizeof(int) * samplesPerCode);

    if (excludeRangeIndex1 < 2)
    {
      for (size_t i = excludeRangeIndex2; i <= (samplesPerCode + excludeRangeIndex1); i++)
      {
        codePhaseRange[newerCodeLen] = i;
        newerCodeLen++;
      }
    }
    else if (excludeRangeIndex2 >= samplesPerCode)
    {
      for (size_t i = (excludeRangeIndex2 - samplesPerCode); i <= excludeRangeIndex1; i++)
      {
        codePhaseRange[newerCodeLen] = i;
        newerCodeLen++;
      }
    }
    else
    {
      int ix = 0;
      for (size_t i = 0; i <= excludeRangeIndex1; i++)
      {
        codePhaseRange[newerCodeLen] = i;
        newerCodeLen++;
      }
      for (size_t i = excludeRangeIndex2; i < samplesPerCode; i++)
      {
        codePhaseRange[newerCodeLen] = i;
        newerCodeLen++;
      }
    }

    // Find the second highest correlation peak in the same freq.bin-- -
    double secondMax = 0.0;
    for (size_t i = 0; i < newerCodeLen; i++)
    {
      if (secondMax < results[frequencyBinIndex][codePhaseRange[i]])
      {
        secondMax = results[frequencyBinIndex][codePhaseRange[i]];
      }
    }
    delete[] codePhaseRange;

    double secondPeakSize = secondMax;

    // -- - Store result---------------------------------------------------- -
    acqResults->peakMetric[PRNIdx] = peakSize / secondPeakSize;

    // If the result is above threshold, then there is a signal ...
    if ((peakSize / secondPeakSize) > settings->acqThreshold)
    {
      *sigFound = true;
      //// Fine resolution frequency search ====================================== =

      // -- - Generate 10msec long C / A codes sequence for given PRN--------
      double* caCode = new double[1023];
      generateCAcode(settings, caCode, PRNIdx);

      int count = 1;
      for (size_t i = 0; i < 10 * samplesPerCode; i++)
      {
        int codeValueIndex = (int)floor((tS * count) / (1.0 / settings->codeFreqBasis));
        codeValueIndex = (int)fmod(codeValueIndex, 1023);
        double longCaCode = caCode[codeValueIndex];

        // (Using detected C / A code phase)
        xCarrier[i] = signal0DC[i + codePhase] * longCaCode;
        count++;
      }
      delete[] caCode;

      // -- - Compute the magnitude of the FFT, find maximumand the
      // associated carrier frequency
      // -- - Find the next highest power of two and increase by 8x--------
      memset(xCarrierIn, 0, sizeof(fftw_complex) * fftNumPts);
      for (size_t i = 0; i < 10 * samplesPerCode; i++)
      {
        xCarrierIn[i][0] = xCarrier[i];
      }
      fft(p4);

      int fftMaxIndex = 0;
      double fftMaxValue = 0.0;
      for (size_t i = 5 - 1; i < (uniqFftPts - 5); i++)
      {
        double absFftXCarrier = xCarrierOut[i][0] * xCarrierOut[i][0] + xCarrierOut[i][1] * xCarrierOut[i][1];
        if (absFftXCarrier > fftMaxValue)
        {
          fftMaxValue = absFftXCarrier;
          fftMaxIndex = i - 4;
        }
      }

      for (size_t i = 0; i < uniqFftPts; i++)
      {
        fftFreqBins[i] = i * settings->samplingFreq / fftNumPts;
      }

      // -- - Save properties of the detected satellite signal------------ -
      acqResults->carrFreq[PRNIdx] = fftFreqBins[fftMaxIndex];
      acqResults->codePhase[PRNIdx] = codePhase;
    }
    else
    {
      // -- - No signal with this PRN--------------------------------------
      //disp(".");
    }// if (peakSize / secondPeakSize) > settings->acqThreshold
  // for PRN = satelliteList

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_destroy_plan(p4);

  delete[] caCodeFreqDomIn;
  delete[] caCodeFreqDom;
  delete[] acqRes1;
  delete[] acqRes2;
  delete[] acqRes1Temp;
  delete[] acqRes2Temp;
  delete[] xCarrierOut;
  delete[] xCarrierIn;
  delete[] fftFreqBins;
  delete[] xCarrier;

  //=== Acquisition is over ==================================================

  for (size_t i = 0; i < numberOfFrqBins; i++)
  {
    delete[] results[i];
  }
  delete[] results;

  for (size_t m = 0; m < (size_t)floor(settings->acqMS/2.0); m++)
  {
    for (size_t i = 0; i < numberOfFrqBins; i++)
    {
      delete[] IQfreqDom1In[m][i];
      delete[] IQfreqDom2In[m][i];
      delete[] IQfreqDom1  [m][i];
      delete[] IQfreqDom2  [m][i];
    }

    delete[] IQfreqDom1In[m];
    delete[] IQfreqDom2In[m];
    delete[] IQfreqDom1  [m];
    delete[] IQfreqDom2  [m];
  }
  delete[] IQfreqDom1In;
  delete[] IQfreqDom2In;
  delete[] IQfreqDom1;
  delete[] IQfreqDom2;
  delete[] IQfreqDom1Loc;
  delete[] IQfreqDom2Loc;

}

 /***********************************
 * Public Function Definitions
 ***********************************/

void modulation(
  double* longSignal,
  const Settings_t* settings
)
{
  double tS = 1.0 / settings->samplingFreq;
  for (size_t i = 0; i < ((size_t)BUFFLENGTH >> 1); i++)
  {
    double cosCarr = cos(2*M_PI*settings->iF*i*tS);
    double sinCarr = sin(2*M_PI*settings->iF*i*tS);
    double IFDblData = cosCarr * longSignal[2 * i + 0] - sinCarr * longSignal[2 * i + 1];
    longSignal[i] = IFDblData;
  }
}

bool acquisition(
  const double* longSignal,
  const Settings_t* settings,
  AcqResults_t* acqResults
)
{
#define FUNC_PREFIX             "acquisition(): "

  //
    // Function performs cold start acquisition on the collected "data".It
    // searches for GPS signals of all satellites, which are listed in field
    // "acqSatelliteList" in the settings structure.Function saves code phase
    // andfrequency of the detected signals in the "acqResults" structure.
    //
    //acqResults = acquisition(longSignal, settings)
    //
    //Inputs:
  // longSignal - 11 ms of raw signal from the front - end
    // settings - Receiver settings->Provides information about
    // sampling and intermediate frequencies and other
    // parameters including the list of the satellites to
    // be acquired.
    // Outputs:
  // acqResults - Function saves code phases and frequencies of the
    // detected signals in the "acqResults" structure.The
    // field "carrFreq" is set to 0 if the signal is not
    // detected for the given PRN number.

    // Find number of samples per spreading code

  bool sigFound = false;
  bool sigFoundAny[32] = { false };

  /**
  */
  double samplesPerCode = round(settings->samplingFreq /
    (settings->codeFreqBasis / settings->codeLength));

  // Create two 1msec vectors of data to correlate with and one with zero DC
  double* signal0DC = new double [BUFFLENGTH];
  subtractMean(signal0DC, longSignal, BUFFLENGTH>>1);

  // Find sampling period
  double tS = 1.0 / settings->samplingFreq ;
  double tC = 1.0 / settings->codeFreqBasis;

  // Generate all C / A codes and sample them according to the sampling freq.
  double** caCodeTables = new double* [32];
  for (int i = 0; i < 32; ++i)
  {
    caCodeTables[i] = new double[samplesPerCode];
  }
  makeCaTable(settings, caCodeTables, 32, samplesPerCode, tS, tC);


  /** Call acquisition function here in multiple threads
  */

  fftw_make_planner_thread_safe();

  // Allocate memory for threads
  thread parallelAcquisitionThreads[32];

  // possible values: 32,16,8,4,2,1
  // reduce 8 to 1 to disable multi threading
  int numParallelThreads = 4;
  int totalThreads = 32 / numParallelThreads;
  for (size_t i = 0; i < totalThreads; i++)
  {
    for (size_t j = 0; j < numParallelThreads; j++)
    {
      int PRNIdx = i * numParallelThreads + j;

      parallelAcquisitionThreads[PRNIdx] = thread(acquisitionCore,
        settings,
        acqResults,
        PRNIdx,
        signal0DC,
        longSignal,
        caCodeTables,
        &sigFoundAny[i]
      );
    }

    for (size_t j = 0; j < numParallelThreads; j++)
    {
      int PRNIdx = i * numParallelThreads + j;
      parallelAcquisitionThreads[PRNIdx].join();
    }
  }

  for (size_t i = 0; i < 32; i++)
  {
    sigFound |= sigFoundAny[i];
  }
  

  for (size_t i = 0; i < 32; i++)
  {
    delete[] caCodeTables[i];
  }
  delete[] caCodeTables;
  delete[] signal0DC;

  return sigFound;
}
