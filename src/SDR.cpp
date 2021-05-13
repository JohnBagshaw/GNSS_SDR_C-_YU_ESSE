/**
 * @file SDR.cpp
 *
 * @brief Main source file for SDR functions in C++
 *
 * Project Title: GNSS-R SDR
 * Author: John Bagshaw
 * Co-author: Surabhi Guruprasad
 * Contact: jotshaw@yorku.ca
 * Supervisor: Prof. Sunil Bisnath
 * Project Manager: Junchan Lee
 * Institution: Lassonde School of Engineering, York University, Canada.
 **/

#include <stdio.h>
#include <iostream>
#include "settings.hpp"
#include "tracking.hpp"
#include "acquisition.hpp"
#include "half_bit_acquisition.hpp"
#include "result_read_write.hpp"
#include <thread>
#include "weak_acquisition.hpp"
#include "ddm_processing.h"

#define _USE_MATH_DEFINES
#include "structures.hpp"
#include "timing.hpp"

#define _CRT_SECURE_NO_WARNINGS

using namespace std;
/***********************************
* Defines
***********************************/

/***********************************
* Macros
***********************************/

/***********************************
* Type Declarations
***********************************/

/***********************************
* Static Member Definitions
***********************************/
size_t mainFunctionExecutionTime;
size_t demodulationFunctionExecutionTime;
size_t normalAcquisitionFunctionExecutionTime;
size_t weakAcquisitionFunctionExecutionTime;
size_t halfBitAcquisitionFunctionExecutionTime;
size_t trackingFunctionExecutionTime;

// Static memory which is being used by tracking part.
static double s_absoluteSample [NUM_CHANNELS][37000];
static double s_codeFreq       [NUM_CHANNELS][37000];
static double s_carrFreq       [NUM_CHANNELS][37000];
static double s_I_P            [NUM_CHANNELS][37000];
static double s_I_E            [NUM_CHANNELS][37000];
static double s_I_L            [NUM_CHANNELS][37000];
static double s_Q_E            [NUM_CHANNELS][37000];
static double s_Q_P            [NUM_CHANNELS][37000];
static double s_Q_L            [NUM_CHANNELS][37000];
static double s_dllDiscr       [NUM_CHANNELS][37000];
static double s_dllDiscrFilt   [NUM_CHANNELS][37000];
static double s_pllDiscr       [NUM_CHANNELS][37000];
static double s_pllDiscrFilt   [NUM_CHANNELS][37000];


static void postProcessing(
  const Settings_t* settings
);



/***********************************
* Entry Point Function
***********************************/
int main()
{

#define FUNC_PREFIX "main(): "
  Tic(1);

  /** reset the timer counters...
  *
  */
  mainFunctionExecutionTime               = 0;
  demodulationFunctionExecutionTime       = 0;
  normalAcquisitionFunctionExecutionTime  = 0;
  weakAcquisitionFunctionExecutionTime    = 0;
  halfBitAcquisitionFunctionExecutionTime = 0;
  trackingFunctionExecutionTime           = 0;

  /** Define main strcutures here which are being used in this application
  *
  *
  * 
  *
  */

  Settings_t globalSettings;
  int startFileNumber = 1;
  int endFileNumber   = 1;

  const int fileNameLen = 200;
  char fileName[fileNameLen];

  for (int fileNumber = startFileNumber; fileNumber <= endFileNumber; fileNumber++)
  {

    /** Populate input binary file name to the variable
    * this file contains raw data which is processed by this application
    *
    */

    memset(fileName, 0, sizeof(char) * fileNameLen);
    //sprintf(fileName, "%sWoodbine_cf7_%d.bin", ".\\in_data\\", fileNumber);
    sprintf(fileName, "%ssample_nottochange_ch1_fileN_%d.bin", ".\\in_data\\", fileNumber);
    printf("Processing File %s \n", fileName);



    /** Initialize settings
    *
    */

    initSettings(fileName, &globalSettings);


    /** Call main post processing function and calculate time spent in it
    *
    */

    if (globalSettings.gnssStart == 1)
    {
      postProcessing(&globalSettings);
    }
    else if (globalSettings.gnssStart == 2)
    {
      DdmProcessing(&globalSettings);
      DelayDopplerMap(&globalSettings);
    }
    else
    {
      disp("No proper gnssStart option specified");
    }
  }


  Toc(1);
  mainFunctionExecutionTime = Delta(1);

  /** Print different execution time stats collected during program execution
  *
  */
  printf("Total Run Time:                 %6.5lf seconds\n", (double)mainFunctionExecutionTime * 1.0e-6);
  printf("Demodulation Run Time:          %6.5lf seconds\n", (double)demodulationFunctionExecutionTime * 1.0e-6);

  switch (globalSettings.acqMethod)
  {
  case SignalAcquisitionMethod_t::NORMAL_SIGNAL_ACQUISITION:
  {
    printf("Normal Acquisition Run Time:    %6.5lf seconds\n", (double)normalAcquisitionFunctionExecutionTime * 1.0e-6);
    break;
  }
  case SignalAcquisitionMethod_t::WEAK_SIGNAL_ACQUISITION:
  {
    printf("Weak Acquisition Run Time:      %6.5lf seconds\n", (double)weakAcquisitionFunctionExecutionTime * 1.0e-6);
    break;
  }
  case SignalAcquisitionMethod_t::HALF_BIT_METHOD:
  {
    printf("Half Bit Acquisition Run Time:  %6.5lf seconds\n", (double)halfBitAcquisitionFunctionExecutionTime * 1.0e-6);
    break;
  }
  default:
    break;
  }

  printf("Tracking Run Time:              %6.5lf seconds\n", (double)trackingFunctionExecutionTime * 1.0e-6);

  return 0;
}

/***********************************
* Static Function Definitions
***********************************/
void postProcessing(
  const Settings_t* settings
) 
{
#define FUNC_PREFIX             "postProcessing(): "

    //// Initialization ======================================================== =
  disp("Starting postProcessing()...");


  /** Define data structures here...
  *
  */
  double snr = 0.0;

  WeakAcqSettings_t weakAcqSettings;

  Channel_t inChannels[NUM_CHANNELS];
  Channel_t outChannels[NUM_CHANNELS];

  AcqResults_t acqResults;
  TrackResults_t trackResults[NUM_CHANNELS];
  IsMatEqualAcqResults_t   isMatEqualAcqResults;
  IsMatEqualTrackResults_t isMatEqualTrackResults[NUM_CHANNELS];

  

  /** Initialize some pointers here
  *
  */

  for (int chIdx = 0; chIdx < NUM_CHANNELS; chIdx++)
  {
    trackResults[chIdx].absoluteSample = s_absoluteSample[chIdx];
    trackResults[chIdx].codeFreq       = s_codeFreq      [chIdx];
    trackResults[chIdx].carrFreq       = s_carrFreq      [chIdx];
    trackResults[chIdx].I_P            = s_I_P           [chIdx];
    trackResults[chIdx].I_E            = s_I_E           [chIdx];
    trackResults[chIdx].I_L            = s_I_L           [chIdx];
    trackResults[chIdx].Q_E            = s_Q_E           [chIdx];
    trackResults[chIdx].Q_P            = s_Q_P           [chIdx];
    trackResults[chIdx].Q_L            = s_Q_L           [chIdx];
    trackResults[chIdx].dllDiscr       = s_dllDiscr      [chIdx];
    trackResults[chIdx].dllDiscrFilt   = s_dllDiscrFilt  [chIdx];
    trackResults[chIdx].pllDiscr       = s_pllDiscr      [chIdx];
    trackResults[chIdx].pllDiscrFilt   = s_pllDiscrFilt  [chIdx];
  }



  /** Open raw data file here so that processing can be done on samples stored in it
  *
  */

  FILE *fid = fopen(settings->fileName, "rb");
  if (fid == NULL)
  {
    disp("File Could not be opened");
  }
  else
  {


    /** File opened successfully, now data can be acquired...
    *
    *
    */

    bool acqResultValid = false;



    /** Do the acquisition here
    *
    *
    */
    if (settings->skipAcquisition == 0)
    {

      /** Define data buffers and read input from the file.
      * Convert int16 data to double data here
      *
      */

      // Find number of samples per spreading 
      double samplesPerCode = round(settings->samplingFreq / (settings->codeFreqBasis / settings->codeLength));
      int16_t *data    = new int16_t [(int)BUFFLENGTH]; // 420112
      double *dataDbl = new double   [(int)BUFFLENGTH]; // 420112


      // Move the starting point of processing. Can be used to start the
      // signal processing at any point in the data record(e.g.good for long
      // records or for signal processing in blocks).
      // Read data for acquisition. 11ms of signal are needed for the fine frequency 
      fseek(fid, 0, SEEK_SET);
      size_t readFlag = fread(data, sizeof(int16_t), (unsigned int)BUFFLENGTH, fid);
      for (size_t i = 0; i < BUFFLENGTH; i++)
      {
        dataDbl[i] = data[i];
      }
      free(data);



      /** Call demodulation function here
      *
      *
      */
      Tic(2);
      disp("Doing the demodulation...");
      modulation(dataDbl, settings);
      Toc(2);
      demodulationFunctionExecutionTime = Delta(2);


      // Signal power
      double Ps = 0.0;
      for (size_t i = 0; i < (size_t)(BUFFLENGTH >> 1); i++)
      {
        Ps += (dataDbl[i] * dataDbl[i]);
      }
      double timeInSeconds = ((BUFFLENGTH >> 1) / settings->samplingFreq);
      Ps /= timeInSeconds;
     
      /* Calculate noise power */
      // Pnoise = kTB = 1.38e–23 * 290 * 2e6
      double No = 8.004e-15;

      /* Calculate SNR now*/
      snr = 10.0 * log10(Ps / No);

      

      /** Call acquisition function here
      *
      *
      */

      switch (settings->acqMethod)
      {


      case SignalAcquisitionMethod_t::NORMAL_SIGNAL_ACQUISITION:
      {

        /** Call normal acquisition function here
        *
        *
        */

        Tic(3);
        disp("Acquiring satellites...");
        acqResultValid = acquisition(
          dataDbl, 
          settings, 
          &acqResults
        );
        Toc(3);
        normalAcquisitionFunctionExecutionTime = Delta(3);

        break;
      }


      case SignalAcquisitionMethod_t::WEAK_SIGNAL_ACQUISITION:
      {

        /** Call normal acquisition function here
        *
        *
        */

        Tic(4);
        disp("Acquiring satellites through weak acquisition...");
        weakAcqSettings.PIT = 10;
        weakAcqSettings.totalSamples = readFlag >> 1;
        fftw_complex* baseBandData = new fftw_complex[weakAcqSettings.totalSamples];
        double tS = 1.0 / settings->samplingFreq;
        for (size_t i = 0; i < ((size_t)BUFFLENGTH >> 1); i++)
        {
          double cosCarr = cos(-2.0 * M_PI * settings->iF * i * tS);
          double sinCarr = sin(-2.0 * M_PI * settings->iF * i * tS);

          baseBandData[i][0] = cosCarr * dataDbl[i];
          baseBandData[i][1] = sinCarr * dataDbl[i];
        }


        acqResultValid =
          weak_acquisition(
            baseBandData,
            settings,
            &weakAcqSettings,
            (WeakAcqResults_t*)&acqResults
          );

        delete[] baseBandData;
        Toc(4);
        weakAcquisitionFunctionExecutionTime = Delta(4);

        break;
      }

      case SignalAcquisitionMethod_t::HALF_BIT_METHOD:
      {
        disp("Acquiring satellites...");
        Tic(5);
        acqResultValid = half_bit_acquisition(
          dataDbl,
          settings,
          &acqResults,
          (int)(BUFFLENGTH >> 1)
        );
        Toc(5);
        halfBitAcquisitionFunctionExecutionTime = Delta(5);

        break;
      }

      default: 
        cout << "No valid choice made for acquisition method" << endl;
        break;

      }


    }




    /** Initialize channelsand prepare for the run
    *
    *
    */

    if (acqResultValid)
    {

      // Start further processing only if a GNSS signal was acquired(the
      // field FREQUENCY will be set to 0 for all not acquired signals)

      preRun(
        &acqResults,
        settings, 
        inChannels
      );

      showChannelStatus(
        settings, 
        inChannels
      );

    }
    else
    {

      // No satellites to track, exit
      disp("No GNSS signals detected, signal processing finished.");
      return;

    }



    /** Track the signal
    *
    *
    */
    Tic(6);
    // Process all channels for given data block
    tracking(
      fid, 
      inChannels, 
      settings, 
      trackResults, 
      outChannels
    );
   Toc(6);
   trackingFunctionExecutionTime = Delta(6);
   disp("Post processing of the signal is over.");



   /** Save Acqusition and tracking results to a MAT file for comparison to be done in MATLAB
   *
   *
   */


   // Auto save the acquisition& tracking results to a file to allow
   // running the positioning solution afterwards.
   disp("Saving Acq & Tracking results to file trackingResults_cpp.mat");
   saveResults(
     &acqResults, 
     &trackResults[0], 
     settings,
     false
   );


   /** Read results from MATLAB and draw comparison with results we got from here...
   *
   *
   */
   /*compareAcquisitionResults(
     settings, 
     &acqResults, 
     &isMatEqualAcqResults
   );
   compareTrackingResults(
     &trackResults[0], 
     &isMatEqualTrackResults[0], 
     settings
   );*/


   /* Display signal to noise ratio here */
   disp("SNR (dBW) : %4.3f\n", snr);

   // Close the file as processing has finished...
   fclose(fid);
  }
}
