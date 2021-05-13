/**
 * @file result_read_write.cpp
 *
 * @brief source file to read/write SDR function results in C++
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
#include "result_read_write.hpp"
#include "mat.h"
#include "matcreat.hpp"
#include <math.h>
#include <algorithm>

 /***********************************
 * Defines
 ***********************************/

 /***********************************
 * Macros
 ***********************************/
#define compareArray(arrayName, arrStr)   { for (size_t chIdx = 0; chIdx < NUM_CHANNELS; chIdx++) \
                                            { \
                                              for (size_t idx = 0; idx < settings->msToProcess; idx++) \
                                              { \
                                                arr[chIdx * NUM_CHANNELS + idx] = trackResults[chIdx].arrayName[idx]; \
                                              } \
                                            } \
                                            writeDoubleValueToMatFile(pmat, &arr[0], NUM_CHANNELS, settings->msToProcess, false, arrStr); \
                                          }\

 /***********************************
 * Static Variables
 ***********************************/

 /***********************************
 * Static Function Definitions
 ***********************************/

static bool compareDoubleArray(
  const double* arr1,
  const double* arr2,
  size_t len
)
{
  bool isSimilar = false;
  double error = 0.0;

  for (size_t i = 0; i < len; i++)
  {
    error += std::abs(arr1[i] - arr2[i]);
  }

  if (error < (double)1.0e-5)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/***********************************
* Public Function Definitions
***********************************/
void showChannelStatus(const Settings_t* settings, const Channel_t* channel)
{
#define FUNC_PREFIX             "showChannelStatus(): "

  int channelNr = 0;
  disp("\n*=========*=====*===============*===========*=============*========*\n");
  disp("| Channel  | PRN  | Frequency  |  Doppler   | Code Offset | Status  |\n");
  disp("*=========*=====*===============*===========*=============*========*\n");

  for (channelNr = 1; channelNr <= settings->numberOfChannels; channelNr++)
  {
#define FUNC_PREFIX "ChannelStatus(): "
    if (channel[channelNr - 1].status != '-')
    {
      disp("|      %2d | %3.0f |  %2.2lf    |   %5.0lf   |    %6lf     |     %c  |\n",
        channelNr,
        channel[channelNr - 1].PRN + 1,
        channel[channelNr - 1].acquiredFreq,
        channel[channelNr - 1].acquiredFreq - settings->iF,
        channel[channelNr - 1].codePhase,
        channel[channelNr - 1].status);
    }
    else
    {
      disp("|      %2d | --- |  ------------ |   -----   |    ------   |   Off  |\n",
        channelNr);
    }
  }
}
void preRun(
  AcqResults_t* acqResults,
  const Settings_t* settings,
  Channel_t* channel
)
{
  int i, j;
  int chIdx;
  int n = 32;
  int minChannels;
  int validChannels;
  int sortedPRNnUmbers[32] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31 };
  /*  PRN number of the tracked satellite */
/*  Used as the center frequency of the NCO */
/*  Position of the C/A  start */
/*  Mode/status of the tracking channel */
/*  "-" - "off" - no signal to track */
/*  "T" - Tracking state */
  static const Channel_t r0 = { 0.0,   /* PRN */
    0.0,                               /* acquiredFreq */
    0.0,                               /* codePhase */
    '-'                                /* status */
  };

  // Initialize to default values
  for (chIdx = 0; chIdx < settings->numberOfChannels; chIdx++) {
    channel[chIdx] = r0;
  }

  double* localPeakMetric = new double[32];
  for (i = 0; i < 32; ++i)
  {

    localPeakMetric[i] = acqResults->peakMetric[i];
  }
  // Sort peaks to find strongest signals, keep the peak index information
  for (i = 0; i < 32; ++i)
  {
    for (j = i + 1; j < 32; ++j)
    {

      if (localPeakMetric[i] < localPeakMetric[j])
      {
        double a = localPeakMetric[i];
        localPeakMetric[i] = localPeakMetric[j];
        localPeakMetric[j] = a;

        a = sortedPRNnUmbers[i];
        sortedPRNnUmbers[i] = sortedPRNnUmbers[j];
        sortedPRNnUmbers[j] = a;
      }
    }
  }

  free(localPeakMetric);

  /* --- Load information about each satellite -------------------------------- */
  /*  Maximum number of initialized channels is number of detected signals, but */
  /*  not more as the number of channels specified in the settings. */

  validChannels = 0;
  for (chIdx = 0; chIdx < 32; chIdx++)
  {
    if (acqResults->carrFreq[chIdx])
    {
      validChannels++;
    }

    if (acqResults->codePhase[chIdx])
    {
      acqResults->codePhase[chIdx] += 1;
    }
  }

  minChannels = std::min(settings->numberOfChannels, validChannels);

  for (chIdx = 0; chIdx < minChannels; chIdx++)
  {
    channel[chIdx].PRN          = sortedPRNnUmbers[chIdx];
    channel[chIdx].acquiredFreq = acqResults->carrFreq[sortedPRNnUmbers[chIdx]];
    channel[chIdx].codePhase    = acqResults->codePhase[sortedPRNnUmbers[chIdx]];
    channel[chIdx].status       = 'T';
  }
}

void compareAcquisitionResults(
  const Settings_t* settings,
  AcqResults_t* acqResults,
  IsMatEqualAcqResults_t* isMatEqualAcqResults
)
{
#define FUNC_PREFIX             "compareAcquisitionResults(): "

  disp("=================> MATLAB vs. C++ Comparison Started\n");

  // Open MAT file to read from
  MATFile* pmat;
  printf("Reading file %s...\n\n", settings->readFileName);
  pmat = matOpen(settings->readFileName, "r");

  if (pmat == NULL) {
    printf("Error creating file %s\n", settings->readFileName);
    printf("(Do you have write permission in this directory?)\n");
  }

  // Allocate memory for holding matlab acquisition array results
  double* carrFreqMat = new double[32];
  double* codePhaseMat = new double[32];
  double* peakMetricMat = new double[32];

  // Read acquisition results from MAT file
  readStructureFieldArrayDouble(pmat, "acqResults", "carrFreq", carrFreqMat, 0);
  readStructureFieldArrayDouble(pmat, "acqResults", "codePhase", codePhaseMat, 0);
  readStructureFieldArrayDouble(pmat, "acqResults", "peakMetric", peakMetricMat, 0);

  // Compare...
  isMatEqualAcqResults->isMatEqual_carrFreq = compareDoubleArray(carrFreqMat, acqResults->carrFreq, 32);
  isMatEqualAcqResults->isMatEqual_codePhase = compareDoubleArray(codePhaseMat, acqResults->codePhase, 32);
  isMatEqualAcqResults->isMatEqual_peakMetric = compareDoubleArray(peakMetricMat, acqResults->peakMetric, 32);

  disp("carrFreq  : %s\n", isMatEqualAcqResults->isMatEqual_carrFreq == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
  disp("codePhase : %s\n", isMatEqualAcqResults->isMatEqual_codePhase == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
  disp("peakMetric: %s\n", isMatEqualAcqResults->isMatEqual_peakMetric == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");

  // Free allocated memory
  delete[] carrFreqMat;
  delete[] codePhaseMat;
  delete[] peakMetricMat;


  // Close opened mat file
  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n", settings->readFileName);
  }
  disp("=================> MATLAB vs. C++ Comparison Ended\n");
}


void compareTrackingResults(
  TrackResults_t* trackResults,
  IsMatEqualTrackResults_t* isMatEqualTrackResults,
  const Settings_t* settings
)
{
#define FUNC_PREFIX             "compareTrackingResults(): "

  disp("=================> MATLAB vs. C++ Comparison Started\n");

  // Open MAT file to read from
  MATFile* pmat;
  printf("Creating file %s...\n\n", settings->readFileName);
  pmat = matOpen(settings->readFileName, "r");

  if (pmat == NULL) {
    printf("Error creating file %s\n", settings->readFileName);
    printf("(Do you have write permission in this directory?)\n");
  }

  isMatEqualTrackResults->isMatEqual_status = true;
  isMatEqualTrackResults->isMatEqual_PRN = true;
  isMatEqualTrackResults->isMatEqual_absoluteSample = true;
  isMatEqualTrackResults->isMatEqual_codeFreq = true;
  isMatEqualTrackResults->isMatEqual_carrFreq = true;
  isMatEqualTrackResults->isMatEqual_I_P = true;
  isMatEqualTrackResults->isMatEqual_I_E = true;
  isMatEqualTrackResults->isMatEqual_I_L = true;
  isMatEqualTrackResults->isMatEqual_Q_E = true;
  isMatEqualTrackResults->isMatEqual_Q_P = true;
  isMatEqualTrackResults->isMatEqual_Q_L = true;
  isMatEqualTrackResults->isMatEqual_dllDiscr = true;
  isMatEqualTrackResults->isMatEqual_dllDiscrFilt = true;
  isMatEqualTrackResults->isMatEqual_pllDiscr = true;
  isMatEqualTrackResults->isMatEqual_pllDiscrFilt = true;

  const char *passStatus = "T";
  for (int chIdx = 0; chIdx < NUM_CHANNELS; chIdx++)
  {
    if (trackResults[chIdx].status == *passStatus)
    {
      // Allocate memory for holding matlab acquisition array results
      char* status_mat = new char[1];
      double msToProcess = settings->msToProcess;
      double* PRN_mat = new double[1];
      double* absoluteSample_mat = new double[msToProcess];
      double* codeFreq_mat = new double[msToProcess];
      double* carrFreq_mat = new double[msToProcess];
      double* I_P_mat = new double[msToProcess];
      double* I_E_mat = new double[msToProcess];
      double* I_L_mat = new double[msToProcess];
      double* Q_E_mat = new double[msToProcess];
      double* Q_P_mat = new double[msToProcess];
      double* Q_L_mat = new double[msToProcess];
      double* dllDiscr_mat = new double[msToProcess];
      double* dllDiscrFilt_mat = new double[msToProcess];
      double* pllDiscr_mat = new double[msToProcess];
      double* pllDiscrFilt_mat = new double[msToProcess];

      // Read acquisition results from MAT file
      readStructureFieldArrayINT8  (pmat, "trackResults", "status", (int8_T*)status_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "PRN", PRN_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "absoluteSample", absoluteSample_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "codeFreq", codeFreq_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "carrFreq", carrFreq_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "I_P", I_P_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "I_E", I_E_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "I_L", I_L_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "Q_E", Q_E_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "Q_P", Q_P_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "Q_L", Q_L_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "dllDiscr", dllDiscr_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "dllDiscrFilt", dllDiscrFilt_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "pllDiscr", pllDiscr_mat, chIdx);
      readStructureFieldArrayDouble(pmat, "trackResults", "pllDiscrFilt", pllDiscrFilt_mat, chIdx);


      // Compare...
      isMatEqualTrackResults->isMatEqual_status &= *status_mat == trackResults[chIdx].status;
      isMatEqualTrackResults->isMatEqual_PRN &= *PRN_mat == trackResults[chIdx].PRN;
      isMatEqualTrackResults->isMatEqual_absoluteSample &= compareDoubleArray(absoluteSample_mat, trackResults[chIdx].absoluteSample, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_codeFreq &= compareDoubleArray(codeFreq_mat, trackResults[chIdx].codeFreq, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_carrFreq &= compareDoubleArray(carrFreq_mat, trackResults[chIdx].carrFreq, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_I_P &= compareDoubleArray(I_P_mat, trackResults[chIdx].I_P, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_I_E &= compareDoubleArray(I_E_mat, trackResults[chIdx].I_E, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_I_L &= compareDoubleArray(I_L_mat, trackResults[chIdx].I_L, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_Q_E &= compareDoubleArray(Q_E_mat, trackResults[chIdx].Q_E, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_Q_P &= compareDoubleArray(Q_P_mat, trackResults[chIdx].Q_P, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_Q_L &= compareDoubleArray(Q_L_mat, trackResults[chIdx].Q_L, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_dllDiscr &= compareDoubleArray(dllDiscr_mat, trackResults[chIdx].dllDiscr, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_dllDiscrFilt &= compareDoubleArray(dllDiscrFilt_mat, trackResults[chIdx].dllDiscrFilt, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_pllDiscr &= compareDoubleArray(pllDiscr_mat, trackResults[chIdx].pllDiscr, (int)msToProcess);
      isMatEqualTrackResults->isMatEqual_pllDiscrFilt &= compareDoubleArray(pllDiscrFilt_mat, trackResults[chIdx].pllDiscrFilt, (int)msToProcess);


      disp("trackResults[channel number = %d].status          : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_status == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].PRN             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_PRN == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].absoluteSample  : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_absoluteSample == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].codeFreq        : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_codeFreq == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].carrFreq        : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_carrFreq == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].I_P             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_I_P == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].I_E             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_I_E == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].I_L             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_I_L == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].Q_E             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_Q_E == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].Q_P             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_Q_P == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].Q_L             : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_Q_L == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].dllDiscr        : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_dllDiscr == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].dllDiscrFilt    : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_dllDiscrFilt == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].pllDiscr        : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_pllDiscr == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");
      disp("trackResults[channel number = %d].pllDiscrFilt    : %s\n", chIdx, isMatEqualTrackResults->isMatEqual_pllDiscrFilt == true ? "MAT vs C++ comparison result: PASSED" : "MAT vs C++ comparison result: FAILED");


      // Free allocated memory
      delete[] status_mat;
      delete[] PRN_mat;
      delete[] absoluteSample_mat;
      delete[] codeFreq_mat;
      delete[] carrFreq_mat;
      delete[] I_P_mat;
      delete[] I_E_mat;
      delete[] I_L_mat;
      delete[] Q_E_mat;
      delete[] Q_P_mat;
      delete[] Q_L_mat;
      delete[] dllDiscr_mat;
      delete[] dllDiscrFilt_mat;
      delete[] pllDiscr_mat;
      delete[] pllDiscrFilt_mat;

    }
  }

  // Close opened mat file
  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n", settings->readFileName);
  }

  disp("=================> MATLAB vs. C++ Comparison Ended\n");
}


void saveResults(
  AcqResults_t* acqResults,
  TrackResults_t* trackResults,
  const Settings_t* settings,
  bool disableTrackResultSave
)
{
  MATFile* pmat;
  //settings->fileName
  printf("Creating file %s...\n\n", settings->writeFileName);
  pmat = matOpen(settings->writeFileName, "w");

  if (pmat == NULL) {
    printf("Error creating file %s\n", settings->writeFileName);
    printf("(Do you have write permission in this directory?)\n");
  }

  // save acquisition results in the mat file
  writeDoubleValueToMatFile(pmat, acqResults->carrFreq, 32, 1, false, "carrFreq");
  writeDoubleValueToMatFile(pmat, acqResults->codePhase, 32, 1, false, "codePhase");
  writeDoubleValueToMatFile(pmat, acqResults->peakMetric, 32, 1, false, "peakMetric");

  if (disableTrackResultSave == false)
  {
    // save tracking results in the mat file
    writeCharValueToMatFile(pmat, &trackResults->status, 1, "status");
    writeDoubleValueToMatFile(pmat, &trackResults->PRN, 1, 1, false, "PRN");

    double* arr = new double[NUM_CHANNELS * settings->msToProcess];
    compareArray(codeFreq, "codeFreq");
    compareArray(carrFreq, "carrFreqt");
    compareArray(I_P, "I_P");
    compareArray(I_E, "I_E");
    compareArray(I_L, "I_L");
    compareArray(Q_E, "Q_E");
    compareArray(Q_P, "Q_P");
    compareArray(Q_L, "Q_L");
    compareArray(dllDiscr, "dllDiscr");
    compareArray(dllDiscrFilt, "dllDiscrFilt");
    compareArray(pllDiscr, "pllDiscr");
    compareArray(pllDiscrFilt, "pllDiscrFilt");
    delete[] arr;
  }

  disp("Done\n");
  if (matClose(pmat) != 0) {
    disp("Error closing file %s\n", settings->writeFileName);
  }

}