/**
 * @file settings.cpp
 *
 * @brief source file for SDR function settings in C++
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
#include "settings.hpp"
#include <string.h>

 /***********************************
 * Defines
 ***********************************/

 /***********************************
 * Macros
 ***********************************/

 /***********************************
 * Static Variables
 ***********************************/

 /***********************************
 * Static Function Definitions
 ***********************************/

void initSettings(
  const char* Nfiles,
  Settings_t* settings
)
{

  int i0;
  // 1 means run post processing, 2 means run DDM processing, else exit
  settings->gnssStart = 2;

  // which acquisition algorithm to run
  settings->acqMethod = SignalAcquisitionMethod_t::NORMAL_SIGNAL_ACQUISITION;

  // This value must be atleast 2 and greater other wise there will be issues with normal acquisition
  settings->acqMS = 4;   

  //  Processing settings ====================================================
  //  Number of milliseconds to be processed used 36000 + any transients (see
  //  below - in Nav parameters) to ensure nav subframes are provided
  settings->msToProcess = 36.0;


  // [ms]
  //  Number of channels to be used for signal processing
  settings->numberOfChannels = NUM_CHANNELS;

  //  Move the starting point of processing. Can be used to start the signal
  //  processing at any point in the data record (e.g. for long records). fseek
  //  function is used to move the file read point, therefore advance is byte
  //  based only.
  settings->skipNumberOfBytes = 0.0;

  //  Raw signal file name and other parameter ===============================
  //  This is a "default" name of the data file (signal record) to be used in
  //  the post-processing mode
  memcpy(&settings->fileName[0], &Nfiles[0], strlen(&Nfiles[0]) + 1);
  int idx = 0;
  int len = strlen(&Nfiles[0]) + 1;
  int lastIdx = 0;
  while (idx < len)
  {
    if (settings->fileName[idx] == '\\')
    {
      lastIdx = idx;
    }
    idx++;
  }

  lastIdx++;
  idx = 0;
  while (settings->fileName[idx + lastIdx] != '.')
  {
    settings->readFileName[idx]  = settings->fileName[idx + lastIdx];
    settings->writeFileName[idx] = settings->fileName[idx + lastIdx];
    settings->ddmWriteFileName[idx] = settings->fileName[idx + lastIdx];    
    idx++;
  }

  const char* ddmPrefix = "_ddm_mat.mat";
  const char* matPrefix = "_mat.mat";
  const char* cppPrefix = "_cpp.mat";
  size_t i = idx;
  for (; i < strlen(matPrefix) + idx; i++)
  {
    settings->readFileName[i] = matPrefix[i - idx];
  }
  settings->readFileName[i] = 0;

  i = idx;
  for (; i < strlen(cppPrefix) + idx; i++)
  {
    settings->writeFileName[i] = cppPrefix[i - idx];
  }
  settings->writeFileName[i] = 0;

  i = idx;
  for (; i < strlen(ddmPrefix) + idx; i++)
  {
    settings->ddmWriteFileName[i] = ddmPrefix[i - idx];
  }
  settings->ddmWriteFileName[i] = 0;

  //  Data type used to store one sample
  settings->dataType = DATA_TYPE_INT16;

  settings->iF = 10e6;

  // [Hz]
  settings->samplingFreq = 40e6;

  // [Hz]
  settings->codeFreqBasis = 1.023E+6;

  // [Hz]
  //  Define number of chips in a code period
  settings->codeLength = 1023.0;

  //  Acquisition settings ===================================================
  //  Skips acquisition in the script postProcessing.m if set to 1
  settings->skipAcquisition = 0.0;

  //  List of satellites to look for. Some satellites can be excluded to speed
  //  up acquisition
  for (i0 = 0; i0 < 32; i0++) {
    settings->acqSatelliteList[i0] = 1.0 + (double)i0;
  }

  // [PRN numbers]
  //  Band around IF to search for satellite signal. Depends on max Doppler
  settings->acqSearchBand = 50;

  // [kHz]
  //  Threshold for the signal presence decision rule
  settings->acqThreshold = 2.5; // [MUST NOT CHANGE IT]

  //  Tracking loops settings ================================================
  //  Code tracking loop parameters
  settings->dllDampingRatio = 0.7;
  settings->dllNoiseBandwidth = 2.0;

  // [Hz]
  settings->dllCorrelatorSpacing = 0.5;

  // [chips]
  //  Carrier tracking loop parameters
  settings->pllDampingRatio = 0.7;
  settings->pllNoiseBandwidth = 25.0;

  // [Hz]
  //  Navigation solution settings ===========================================
  //  Period for calculating pseudoranges and position
  settings->navSolPeriod = 500.0;

  // [ms]
  //  Elevation mask to exclude signals from satellites at low elevation
  settings->elevationMask = 10.0;

  // [degrees 0 - 90]
  //  Enable/dissable use of tropospheric correction
  settings->useTropCorr = 1.0;

  //  0 - Off
  //  1 - On
  //  True position of the antenna in UTM system (if known). Otherwise enter
  //  all NaN's and mean position will be used as a reference .
  settings->truePosition.E = -1.0;
  settings->truePosition.N = -1.0;
  settings->truePosition.U = -1.0;

  // [ms] Initial sign. travel time
  settings->startOffset = 68.802;
}
