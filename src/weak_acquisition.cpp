/**
 * @file weak_sig_acquisition.cpp
 *
 * @brief source file for SDR weak signal acquisition algorithm in C++
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
#include <math.h>
#include "gen_code.hpp"
#include <string.h>
#include "constants.h"
#include "result_read_write.hpp"
#include "weak_acquisition.hpp"
#include <thread>
#include <mutex>

 /***********************************
 * Defines
 ***********************************/

 /***********************************
 * Macros
 ***********************************/
std::mutex mtx2;
std::mutex mtx3;

static void fft(
  fftw_plan p
)
{
  

  /*mtx2.lock();*/
  //mtx2.unlock();

  fftw_execute(p);

  //mtx2.lock();
  //mtx2.unlock();
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

void weakAcquisitionCore(
  int                      PRNIdx,
  const fftw_complex      *longSignal,
  const Settings_t        *settings,
  const WeakAcqSettings_t *weakAcqSetting,
        WeakAcqResults_t  * acqResults,
        double ** caCodeTables,
        double samplesPerCode
)
{
#define FUNC_PREFIX "Weak Acquisition Parallel Function(): "


  /** Reset acq results structures
  *
  *
  */
  
  acqResults->carrFreq[PRNIdx] = 0.0;   // Carrier frequencies of detected signals
  acqResults->codePhase[PRNIdx] = 0.0;  // C / A code phases of detected signals
  acqResults->peakMetric[PRNIdx] = 0.0; // Correlation peak ratios of the detected signals

  /** Acquisition parameters
  *
  *
  */


  double numberOfFrqBins = round(settings->acqSearchBand * 2) + 1;
  int Nfd = (int)settings->acqSearchBand / 0.5 + 1;
  int Sblock = (int)floor((samplesPerCode * weakAcqSetting->PIT) / Nfd);
  int Nint = (int)Sblock * Nfd;
  int Nblocks = (int)floor(weakAcqSetting->totalSamples / Sblock) - Nfd;

  /** Reserve memories
  *
  *
  */

  fftw_complex* sigRowFftIn      = new fftw_complex  [Sblock];
  fftw_complex* sigRowFftOut     = new fftw_complex  [Sblock];
  fftw_complex* caCodeFftIn      = new fftw_complex  [Sblock];
  fftw_complex* caCodeFftOut     = new fftw_complex  [Sblock];
  fftw_complex* sigCaCodeFftIn   = new fftw_complex  [Sblock];
  fftw_complex* sigCaCodeFftOut  = new fftw_complex  [Sblock];
  fftw_complex* corrMatrixColIn  = new fftw_complex  [Nfd];
  fftw_complex* corrMatrixColOut = new fftw_complex  [Nfd];
  fftw_complex** corrMatrix      = new fftw_complex *[Nfd];
  fftw_complex** corrFftMatrix   = new fftw_complex *[Nfd];
  fftw_complex** sigMatrix       = new fftw_complex *[Nfd];

        double** caCodeMatrix   = new double* [Nfd];

  for (int i = 0; i < Nfd; ++i)
  {
    corrMatrix    [i] = new fftw_complex [Sblock];
    corrFftMatrix [i] = new fftw_complex [Sblock];
    sigMatrix     [i] = new fftw_complex [Sblock];
    caCodeMatrix  [i] = new double       [Sblock];
  }

  int PRN = settings->acqSatelliteList[PRNIdx] - 1.0;
  int block = 0;
  int blockStart = 0;
  bool foundSignal = false;

  fftw_plan p1 = fftw_plan_dft_1d(Sblock, sigRowFftIn    , sigRowFftOut    , true  ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d(Sblock, caCodeFftIn    , caCodeFftOut    , true  ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p3 = fftw_plan_dft_1d(Sblock, sigCaCodeFftIn , sigCaCodeFftOut , false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p4 = fftw_plan_dft_1d(Nfd,    corrMatrixColIn, corrMatrixColOut, true  ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

  while (block < Nblocks)
  {

    //  2. Define block length and rearrange the signal into the matrix
    int idx = 0;
    for (int row = 0; row < Nfd; row++)
    {
      for (int col = 0; col < Sblock; col++)
      {
        sigMatrix[row][col][0] = longSignal[blockStart + idx][0] + longSignal[blockStart + Sblock + idx][0];
        sigMatrix[row][col][1] = longSignal[blockStart + idx][1] + longSignal[blockStart + Sblock + idx][1];
        idx++;
      }
    }

    //  3. Define local PN code and rearrange it also into the matrix
    idx = 0;
    for (int row = 0; row < Nfd; row++)
    {
      for (int col = 0; col < Sblock; col++)
      {
        caCodeMatrix[row][col] = caCodeTables[PRN][idx];
        idx += 1;
        idx %= (int)samplesPerCode;
      }
    }

    //  4. Do the circular correlation for each of the row
    for (int row = 0; row < Nfd; row++)
    {
      for (int col = 0; col < Sblock; col++)
      {
        sigRowFftIn[col][0] = sigMatrix[row][col][0];
        sigRowFftIn[col][1] = sigMatrix[row][col][1];

        caCodeFftIn[col][0] = caCodeMatrix[row][col];
        caCodeFftIn[col][1] = 0.0;
      }

      fft(p1);
      fft(p2);

      for (int col = 0; col < Sblock; col++)
      {
        caCodeFftOut[col][0] = +1.0 * caCodeFftOut[col][0];
        caCodeFftOut[col][1] = -1.0 * caCodeFftOut[col][1];
      }

      for (int col = 0; col < Sblock; col++)
      {
        sigCaCodeFftIn[col][0] = sigRowFftOut[col][0] * caCodeFftOut[col][0] - sigRowFftOut[col][1] * caCodeFftOut[col][1];
        sigCaCodeFftIn[col][1] = sigRowFftOut[col][0] * caCodeFftOut[col][1] + sigRowFftOut[col][1] * caCodeFftOut[col][0];
      }

      fft(p3);

      for (int col = 0; col < Sblock; col++)
      {
        corrMatrix[row][col][0] = sigCaCodeFftOut[col][0] * (1.0 / Sblock);
        corrMatrix[row][col][1] = sigCaCodeFftOut[col][1] * (1.0 / Sblock);
      }
    }

    //  5. Do the DFT for each column and check if not greater than threshold

    for (int col = 0; col < Sblock; col++)
    {
      for (int row = 0; row < Nfd; row++)
      {
        corrMatrixColIn[row][0] = corrMatrix[row][col][0];
        corrMatrixColIn[row][1] = corrMatrix[row][col][1];
      }

      fft(p4);


      for (int row = 0; row < Nfd; row++)
      {
        corrFftMatrix[row][col][0] = corrMatrixColOut[row][0] * corrMatrixColOut[row][0] +
          corrMatrixColOut[row][1] * corrMatrixColOut[row][1];
      }
    }

    //  6. If not, move one block onto next data and start again
    int freqBin;
    double* colVec = new double[Nfd];
    double maxVal = 0;
    for (int row = 0; row < Nfd; row++)
    {
      maxVal = 0;
      for (int col = 0; col < Sblock; col++)
      {
        if (corrFftMatrix[row][col][0] > maxVal)
        {
          maxVal = corrFftMatrix[row][col][0];
        }
      }

      colVec[row] = maxVal;
    }

    double maxValue = 0.0;
    for (int row = 0; row < Nfd; row++)
    {
      if (colVec[row] > maxValue)
      {
        maxValue = colVec[row];
        freqBin = row;
      }
    }
    delete[] colVec;

    double maxFftAbs = 0, codePhase;

    double* rowVec = new double[Sblock];
    for (int col = 0; col < Sblock; col++)
    {
      maxFftAbs = 0;
      for (int row = 0; row < Nfd; row++)
      {
        if (corrFftMatrix[row][col][0] > maxFftAbs)
        {
          maxFftAbs = corrFftMatrix[row][col][0];
        }
      }

      rowVec[col] = maxFftAbs;
    }

    maxFftAbs = 0;
    for (int col = 0; col < Sblock; col++)
    {
      if (rowVec[col] > maxFftAbs)
      {
        maxFftAbs = rowVec[col];
        codePhase = col;
      }
    }
    delete[] rowVec;

    double samplesPerCodeChip = round(settings->samplingFreq / settings->codeFreqBasis);
    double excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    double excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    double secondMaxFftAbs = 0;
    for (int i = 0; i <= excludeRangeIndex1; i++)
    {
      if (corrFftMatrix[freqBin][i][0] > secondMaxFftAbs)
      {
        secondMaxFftAbs = corrFftMatrix[freqBin][i][0];
      }
    }
    for (int i = excludeRangeIndex2; i < Sblock; i++)
    {
      if (corrFftMatrix[freqBin][i][0] > secondMaxFftAbs)
      {
        secondMaxFftAbs = corrFftMatrix[freqBin][i][0];
      }
    }

    double peakTo2ndPeakRatio = maxFftAbs / secondMaxFftAbs;
    acqResults->peakMetric[PRN] = peakTo2ndPeakRatio;

    //  If the result is above threshold, then there is a signal ...
    if (peakTo2ndPeakRatio > settings->acqThreshold)
    {
      //  Save properties of the detected satellite signal
      acqResults->carrFreq[PRN] = settings->iF + (1 / (weakAcqSetting->PIT * 1e-3 * (freqBin + 1)));
      acqResults->codePhase[PRN] = codePhase;

      foundSignal = true;
      break;
    }

    block = block + 1;
    blockStart = blockStart + Sblock;

  } //  for next block in the matrix

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(p3);
  fftw_destroy_plan(p4);

  if (foundSignal == false)
  {
    mtx3.lock();
    disp(".");
    mtx3.unlock();
  }
  else
  {
    mtx3.lock();
    disp("%d ", PRN);
    mtx3.unlock();
  }

  delete[] sigRowFftIn;
  delete[] sigRowFftOut;
  delete[] caCodeFftIn;
  delete[] caCodeFftOut;
  delete[] sigCaCodeFftIn;
  delete[] sigCaCodeFftOut;
  delete[] corrMatrixColIn;
  delete[] corrMatrixColOut;

  for (int i = 0; i < Nfd; ++i)
  {
    delete[] corrMatrix[i];
    delete[] corrFftMatrix[i];
    delete[] sigMatrix[i];
    delete[] caCodeMatrix[i];
  }
  delete[] corrMatrix;
  delete[] corrFftMatrix;
  delete[] sigMatrix;
  delete[] caCodeMatrix;
}

 /***********************************
 * Public Function Definitions
 ***********************************/

bool weak_acquisition(
  const fftw_complex* longSignal,
  const Settings_t* settings,
  const WeakAcqSettings_t* weakAcqSetting,
  WeakAcqResults_t* acqResults
)
{
#define FUNC_PREFIX             "weak_acquisition(): "

  //  Find number of samples per spreading code
  double samplesPerCode = round(settings->samplingFreq /
    (settings->codeFreqBasis / settings->codeLength));

  //  Find sampling period
  double tS = 1 / settings->samplingFreq;
  double tC = 1.0 / settings->codeFreqBasis;

  //  Generate all C/A codes and sample them according to the sampling freq.
  double** caCodeTables = new double* [32];
  for (int i = 0; i < 32; ++i)
  {
    caCodeTables[i] = new double[samplesPerCode];
  }
  makeCaTable(settings, caCodeTables, 32, samplesPerCode, tS, tC);


  disp("(");

  /** Perform search for all listed PRN numbers
  *
  *
  */  

  fftw_make_planner_thread_safe();

  thread weakAcquisitionParallelThreads[32];
  int numParallelThreads = 32;
  int totalThreads = 32 / numParallelThreads;

  for (size_t i = 0; i < totalThreads; i++)
  {
    for (size_t j = 0; j < numParallelThreads; j++)
    {
      int PRNIdx = i * numParallelThreads + j;

      weakAcquisitionParallelThreads[PRNIdx] = thread(weakAcquisitionCore,
        PRNIdx,
        longSignal,
        settings,
        weakAcqSetting,
        acqResults,
        caCodeTables,
        samplesPerCode
      );
    }

    for (size_t j = 0; j < numParallelThreads; j++)
    {
      int PRNIdx = i * numParallelThreads + j;
      weakAcquisitionParallelThreads[PRNIdx].join();
    }
  }

  // === Acquisition is over ==================================================
  disp(")\n");

  return true;
}
