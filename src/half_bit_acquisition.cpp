/**
 * @file half_bit_acquisition.cpp
 *
 * @brief alternate half-bit acquisition source file for SDR functions in C++
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
#include "half_bit_acquisition.hpp"
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
static std::mutex mtx;

static void fft(
  fftw_plan p
) 
{ 
  fftw_execute(p);
}


 /***********************************
 * Static Variables
 ***********************************/
extern size_t sAcquisitionTime;
extern size_t sAcquisitionOuterLoop;
extern size_t sAcquisitionInnerLoop;
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

static void halfBitAcquisitionCore(
  const Settings_t* settings,
  AcqResults_t* acqResults,
  int    PRNIdx,
  int data_len_ms,
  double* signal1,
  double* signal2,
  double** caCodeTables,
  bool* sigFound
)
{

#define FUNC_PREFIX "Half Bit Acquisition: "
  mtx.lock();
  disp("PRN = %3i\n", PRNIdx + 1);
  mtx.unlock();

  // -- - Initialize arrays to speed up the code------------------------------ -
  // Search results of all frequency bins and code shifts(for one satellite)

  // Number of the frequency bins for the given acquisition band(500Hz steps)
  double samplesPerCode = round(settings->samplingFreq / (settings->codeFreqBasis / settings->codeLength));
  int numberOfFrqBins = round(settings->acqSearchBand * 2) + 1;
  double tS = 1 / settings->samplingFreq;
  double tC = 1.0 / settings->codeFreqBasis;
  int fftNumPts = (int)(8 * pow(2, ceil(log((data_len_ms-1) * samplesPerCode) / log(2))));
  int uniqFftPts = (int)ceil((fftNumPts + 1) / 2.0);

  double* xCarrier = new double[10 * samplesPerCode];
  double* fftFreqBins = new double[uniqFftPts];

  double** resultsEven = new double* [numberOfFrqBins];
  double** resultsOdd  = new double* [numberOfFrqBins];
  double** results     = new double* [numberOfFrqBins];

  for (int i = 0; i < numberOfFrqBins; ++i)
  {
    resultsEven [i] = new double [samplesPerCode];
    resultsOdd  [i] = new double [samplesPerCode];
    results     [i] = new double [samplesPerCode];
  }

  double* signal0DC = new double[data_len_ms * samplesPerCode];

  fftw_complex* caCodeFreqDomIn = new fftw_complex [samplesPerCode];
  fftw_complex* caCodeFreqDom   = new fftw_complex [samplesPerCode];

  fftw_complex* acqResEven = new fftw_complex [samplesPerCode];
  fftw_complex* acqResOdd  = new fftw_complex [samplesPerCode];

  fftw_complex* xCarrierIn  = new fftw_complex [fftNumPts     ];
  fftw_complex* xCarrierOut = new fftw_complex [fftNumPts     ];

  fftw_complex*** IQfreqDomEvenIn = new fftw_complex ** [data_len_ms];
  fftw_complex*** IQfreqDomEven   = new fftw_complex ** [data_len_ms];
  fftw_complex*** IQfreqDomOddIn  = new fftw_complex ** [data_len_ms];
  fftw_complex*** IQfreqDomOdd    = new fftw_complex ** [data_len_ms];

  for (int i = 0; i < data_len_ms; ++i)
  {
    IQfreqDomEvenIn[i] = new fftw_complex * [numberOfFrqBins];
    IQfreqDomEven  [i] = new fftw_complex * [numberOfFrqBins];
    IQfreqDomOddIn [i] = new fftw_complex * [numberOfFrqBins];
    IQfreqDomOdd   [i] = new fftw_complex * [numberOfFrqBins];

    for (int j = 0; j < numberOfFrqBins; ++j)
    {
      IQfreqDomEvenIn[i][j] = new fftw_complex [samplesPerCode];
      IQfreqDomEven  [i][j] = new fftw_complex [samplesPerCode];
      IQfreqDomOddIn [i][j] = new fftw_complex [samplesPerCode];
      IQfreqDomOdd   [i][j] = new fftw_complex [samplesPerCode];
    }
  }

  fftw_complex* IQfreqDomEvenLoc = new fftw_complex[samplesPerCode];
  fftw_complex* IQfreqDomOddLoc  = new fftw_complex[samplesPerCode];

  double scale = 1.0 / samplesPerCode;
  for (int i = 0; i < numberOfFrqBins; ++i)
  {
    for (int j = 0; j < data_len_ms; j++)
    {
      fftw_plan iqFreqDomEvenPlans = fftw_plan_dft_1d((int)samplesPerCode, &IQfreqDomEvenIn[j][i][0], &IQfreqDomEven[j][i][0], true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
      fftw_plan iqFreqDomOddPlans  = fftw_plan_dft_1d((int)samplesPerCode, &IQfreqDomOddIn [j][i][0], &IQfreqDomOdd [j][i][0], true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

      double frqBins = settings->iF - (settings->acqSearchBand / 2) * 1000 + 0.5e3 * (i);
      for (size_t k = 0; k < samplesPerCode; k++)
      {
        int idx = j * samplesPerCode + k;
        double sinCar = sin(frqBins * k * 2 * Pi * tS);
        double cosCar = cos(frqBins * k * 2 * Pi * tS);

        // Generate carrier wave frequency grid(0.5kHz step)-----------
        // -- - Generate local sine and cosine------------------------------
        // -- - Convert the baseband signal to frequency domain-------------
        IQfreqDomEvenIn[j][i][k][0] = sinCar * signal1[idx];
        IQfreqDomEvenIn[j][i][k][1] = cosCar * signal1[idx];
        IQfreqDomOddIn [j][i][k][0] = sinCar * signal2[idx];
        IQfreqDomOddIn [j][i][k][1] = cosCar * signal2[idx];
      }

      fft(iqFreqDomEvenPlans);
      fft(iqFreqDomOddPlans);

      fftw_destroy_plan(iqFreqDomEvenPlans);
      fftw_destroy_plan(iqFreqDomOddPlans);
    }
  }


  fftw_plan p1 = fftw_plan_dft_1d((int)samplesPerCode, caCodeFreqDomIn , caCodeFreqDom, false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d((int)samplesPerCode, IQfreqDomEvenLoc, acqResEven   , false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p3 = fftw_plan_dft_1d((int)samplesPerCode, IQfreqDomOddLoc , acqResOdd    , false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p4 = fftw_plan_dft_1d(fftNumPts          , xCarrierIn      , xCarrierOut  , true  ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

  // -- - Initialize acqResults------------------------------------------------
  
  // Carrier frequencies of detected signals
  acqResults->carrFreq[PRNIdx] = 0.0;
  // C / A code phases of detected signals
  acqResults->codePhase[PRNIdx] = 0.0;
  // Correlation peak ratios of the detected signals
  acqResults->peakMetric[PRNIdx] = 0.0;

  int PRN = (int)settings->acqSatelliteList[PRNIdx];

  for (int i = 0; i < numberOfFrqBins; i++)
  {
    for (int j = 0; j < samplesPerCode; j++)
    {
      resultsEven[i][j] = 0.0;
      resultsOdd [i][j] = 0.0;
    }
  }


    //// Correlate signals ======================================================
    for (size_t i = 0; i < (size_t)samplesPerCode; i++)
    {
      caCodeFreqDomIn[i][0] = caCodeTables[PRNIdx][i];
      caCodeFreqDomIn[i][1] = 0.0;
    }
    fft(p1);

    // Make the correlation for whole frequency band(for all freq.bins)
    for (size_t frqBinIndex = 0; frqBinIndex < numberOfFrqBins; frqBinIndex++)
    {
      // -- - Multiplication in the frequency domain(correlation in time domain)
      for (int ms = 0; ms < data_len_ms; ms++)
      {

        for (size_t i = 0; i < samplesPerCode; i++)
        {
          IQfreqDomEvenLoc[i][0] = (IQfreqDomEven[ms][frqBinIndex][i][0] * caCodeFreqDom[i][0]) - (IQfreqDomEven[ms][frqBinIndex][i][1] * caCodeFreqDom[i][1]);
          IQfreqDomEvenLoc[i][1] = (IQfreqDomEven[ms][frqBinIndex][i][0] * caCodeFreqDom[i][1]) + (IQfreqDomEven[ms][frqBinIndex][i][1] * caCodeFreqDom[i][0]);

          IQfreqDomOddLoc[i][0] = (IQfreqDomOdd[ms][frqBinIndex][i][0] * caCodeFreqDom[i][0]) - (IQfreqDomOdd[ms][frqBinIndex][i][1] * caCodeFreqDom[i][1]);
          IQfreqDomOddLoc[i][1] = (IQfreqDomOdd[ms][frqBinIndex][i][0] * caCodeFreqDom[i][1]) + (IQfreqDomOdd[ms][frqBinIndex][i][1] * caCodeFreqDom[i][0]);
        }
        fft(p2);
        fft(p3);

        // -- - Perform inverse DFT and store correlation results------------
        double maxAcqRes1 = 0.0;
        double maxAcqRes2 = 0.0;

        for (size_t i = 0; i < samplesPerCode; i++)
        {
          acqResEven[i][0] *= scale;
          acqResOdd[i][0] *= scale;
          acqResEven[i][1] *= scale;
          acqResOdd[i][1] *= scale;

          acqResEven[i][0] = acqResEven[i][0] * acqResEven[i][0] + acqResEven[i][1] * acqResEven[i][1];
          acqResOdd [i][0] = acqResOdd[i][0] * acqResOdd[i][0] + acqResOdd[i][1] * acqResOdd[i][1];

          resultsEven[frqBinIndex][i] += acqResEven[i][0];
          resultsOdd [frqBinIndex][i] += acqResOdd [i][0];
        }
      }
    }

     /*
      Check which msec had the greater power and save that, will
      "blend" 1st and 2nd msec but will correct data bit issues
      compare between even and odd sets
     */
      double evenPeakValue = 0.0;
      double oddPeakValue  = 0.0;
      for (size_t f = 0; f < numberOfFrqBins; f++)
      {
        for (size_t i = 0; i < samplesPerCode; i++)
        {
          if (resultsEven[f][i] > evenPeakValue)
          {
            evenPeakValue = resultsEven[f][i];
          }
          if (resultsOdd[f][i] > oddPeakValue)
          {
            oddPeakValue = resultsOdd[f][i];
          }
        }
      }

      for (size_t f = 0; f < numberOfFrqBins; f++)
      {
        for (size_t i = 0; i < samplesPerCode; i++)
        {
          if (evenPeakValue > oddPeakValue)
          {
            results[f][i] = resultsEven[f][i];
          }
          else
          {
            results[f][i] = resultsOdd[f][i];
          }
        }
      }
      if (evenPeakValue > oddPeakValue)
      {
        subtractMean(signal0DC, signal1, data_len_ms * samplesPerCode );
      }
      else
      {
        subtractMean(signal0DC, signal2,data_len_ms * samplesPerCode);
      }
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
      for (size_t i = 0; i < (data_len_ms-1) * samplesPerCode; i++)
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
      for (size_t i = 0; i < (data_len_ms-1) * samplesPerCode; i++)
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
  delete[] acqResEven;
  delete[] acqResOdd;
  delete[] xCarrierOut;
  delete[] xCarrierIn;
  delete[] fftFreqBins;
  delete[] xCarrier;

  //=== Acquisition is over ==================================================

  for (size_t i = 0; i < numberOfFrqBins; i++)
  {
    delete[] resultsEven[i];
    delete[] resultsOdd[i];
    delete[] results[i];
  }
  delete[] results;
  delete[] resultsEven;
  delete[] resultsOdd;

  for (int i = 0; i < data_len_ms; i++)
  {
    for (int j = 0; j < numberOfFrqBins; j++)
    {
      delete[] IQfreqDomEvenIn[i][j];
      delete[] IQfreqDomOddIn[i][j];
      delete[] IQfreqDomEven[i][j];
      delete[] IQfreqDomOdd[i][j];
    }
    delete[] IQfreqDomEvenIn[i];
    delete[] IQfreqDomOddIn[i];
    delete[] IQfreqDomEven[i];
    delete[] IQfreqDomOdd[i];
  }
  delete[] IQfreqDomEvenIn;
  delete[] IQfreqDomOddIn;
  delete[] IQfreqDomEven;
  delete[] IQfreqDomOdd;
  delete[] IQfreqDomEvenLoc;
  delete[] IQfreqDomOddLoc;
  delete[] signal0DC;
}

 /***********************************
 * Public Function Definitions
 ***********************************/

bool half_bit_acquisition(
  const double* longSignal,
  const Settings_t* settings,
  AcqResults_t* acqResults,
  int dataLength
)
{

#define FUNC_PREFIX             "halfBitAcquisition(): "

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
  int data_len_ms = 10;
  if (floor((double)dataLength / samplesPerCode) > 20.0)
  {
    data_len_ms = 10;
  }
  else if (floor((double)dataLength / samplesPerCode) > 12.0)
  {
    data_len_ms = floor((double)dataLength / samplesPerCode) - 10;
  }
  else
  {
    disp("Error: not enough data is available for half bit acquisition method");
    exit(1);
  }

  double* signal1   = new double [data_len_ms * samplesPerCode];
  double* signal2   = new double [data_len_ms * samplesPerCode];

  memcpy(signal1, &longSignal[0],                      sizeof(double) * data_len_ms * samplesPerCode);
  memcpy(signal2, &longSignal[10*(int)samplesPerCode], sizeof(double) * data_len_ms * samplesPerCode);

  // Find sampling period
  double tS = 1 / settings->samplingFreq;
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

  int numParallelThreads = 4;
  int totalThreads = 32 / numParallelThreads;
  for (size_t i = 0; i < totalThreads; i++)
  {
    for (size_t j = 0; j < numParallelThreads; j++)
    {
      int PRNIdx = i * numParallelThreads + j;

      parallelAcquisitionThreads[PRNIdx] = thread(halfBitAcquisitionCore,
        settings,
        acqResults,
        PRNIdx,
        data_len_ms,
        signal1,
        signal2,
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
  delete[] signal1;
  delete[] signal2;

  return sigFound;
}
