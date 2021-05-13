/**
 * @file tracking.cpp
 *
 * @brief source file for SDR tracking functions in C++
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
#include "tracking.hpp"
#include "result_read_write.hpp"
#include "gen_code.hpp"

#include <math.h>
#include <stdint.h>
#include <string>
#include "timing.hpp"
#include <thread>
#include <mutex>

 /***********************************
 * Defines
 ***********************************/
#define _CRT_DISABLE_PERFCRIT_LOCKS
 /***********************************
 * Macros
 ***********************************/

 /***********************************
 * Static Variables
 ***********************************/
mutex mtx1;
 /***********************************
 * Static Function Definitions
 ***********************************/
void trackingCore(
int channelNr,
const Settings_t *settings,
FILE* fid,
const Channel_t* InChannels,
TrackResults_t* trackResults,
Channel_t* outChannels
)
{
#define FUNC_PREFIX "TrackingCore Parallel Function: "

  int loopCnt = 0;

  //// Initialize tracking variables ==========================================

  double codePeriods = settings->msToProcess;// For GPS one C / A code is one ms

  // -- - DLL variables--------------------------------------------------------
  // Define early - late offset(in chips)
  double earlyLateSpc = settings->dllCorrelatorSpacing;

  // Summation interval
  double PDIcode = 0.001;

  // Calculate filter coefficient values

  double LBW = settings->dllNoiseBandwidth;
  double zeta = settings->dllDampingRatio;
  double k = 1.0;

  // Solve natural frequency
  double Wn = LBW * 8 * zeta / ((4 * zeta * zeta) + 1);

  // solve for t1& t2
  double tau1code = k / (Wn * Wn);
  double tau2code = 2.0 * zeta / Wn;

  // -- - PLL variables--------------------------------------------------------
    // Summation interval
  double PDIcarr = 0.001;

  LBW = settings->pllNoiseBandwidth;
  zeta = settings->pllDampingRatio;
  k = 0.25;

  // Solve natural frequency
  Wn = LBW * 8 * zeta / ((4 * zeta * zeta) + 1);

  // solve for t1& t2
  double tau1carr = k / (Wn * Wn);
  double tau2carr = 2.0 * zeta / Wn;


  // Only process if PRN is non zero(acquisition was successful)
  if ((InChannels[channelNr].PRN + 1) != 0)
  {

    // Save additional information - each channel's tracked PRN
    trackResults[channelNr].PRN = InChannels[channelNr].PRN + 1;

    // Get a vector with the C / A code sampled 1x / chip
    double caCode1[1023];
    double caCode[1025];
    generateCAcode(settings, caCode1, trackResults[channelNr].PRN - 1);

    // Then make it possible to do earlyand late versions
    caCode[0] = caCode1[1022];
    caCode[1024] = caCode1[0];
    for (size_t i = 0; i < 1023; i++)
    {
      caCode[i + 1] = caCode1[i];
    }

    // -- - Perform various initializations------------------------------
      // define initial code frequency basis of NCO
    double codeFreq = settings->codeFreqBasis;

    // define residual code phase(in chips)
    double remCodePhase = 0.0;

    // define carrier frequency which is used over whole tracking period
    double carrFreq = InChannels[channelNr].acquiredFreq;
    double carrFreqBasis = InChannels[channelNr].acquiredFreq;

    // define residual carrier phase
    double remCarrPhase = 0.0;

    // code tracking loop parameters
    double oldCodeNco = 0.0;
    double oldCodeError = 0.0;

    // carrier / Costas loop parameters
    double oldCarrNco = 0.0;
    double oldCarrError = 0.0;

    int fileReadLoc = settings->skipNumberOfBytes + InChannels[channelNr].codePhase;

    //=== Process the number of specified code periods ================
    for (loopCnt = 0; loopCnt < (int)codePeriods; loopCnt++)
    {
      mtx1.lock();
      disp("Channel Number: %i / %i, Loop Counter: %i / %i\n", channelNr + 1, settings->numberOfChannels, loopCnt + 1, (int)codePeriods);
      mtx1.unlock();

      // //Read next block of data------------------------------------------------
      // Find the size of a "block" or code period in whole samples

      // Update the phasestep based on code freq(variable) and
      //sampling frequency(fixed)
      double codePhaseStep = codeFreq / settings->samplingFreq;

      int blksize = ceil((settings->codeLength - remCodePhase) / codePhaseStep);

      // Read in the appropriate number of samples to process this
      // interation

      mtx1.lock();
      // Move the starting point of processing.Can be used to start the
      // signal processing at any point in the data record(e.g. for long
      // records).In addition skip through that data file to start at the
      // appropriate sample(corresponding to code phase).Assumes sample
      // type is schar(or 1 byte per sample)
      if (fileReadLoc == 0)
      {
        fileReadLoc = 1;
      }

      fseek(fid, fileReadLoc-1, SEEK_SET);
      int16_t* rawSignal = (int16_t*)malloc(sizeof(int16_t) * blksize);
      size_t readCount = fread(rawSignal, sizeof(int16_t), (unsigned int)blksize, fid);
      fileReadLoc += (readCount * sizeof(int16_t));
      mtx1.unlock();

      // If did not read in enough samples, then could be out of
        // data - better exit
      if (readCount != blksize)
      {
        disp("Not able to read the specified number of samples  for tracking, exiting!");
        fclose(fid);
        return;
      }

      // //Set up all the code phase tracking information------------------------ -

      // Define index into early code vector
      int tcode;
      double tcodeInit;
      double* earlyCode = (double*)malloc(sizeof(double) * blksize);
      tcodeInit = (remCodePhase - earlyLateSpc);
      for (size_t i = 0; i < blksize; i++)
      {
        tcode = (int)ceil(tcodeInit + i * codePhaseStep);
        earlyCode[i] = caCode[tcode];
      }

      double* lateCode = (double*)malloc(sizeof(double) * blksize);
      tcodeInit = (remCodePhase + earlyLateSpc);
      for (size_t i = 0; i < blksize; i++)
      {
        tcode = (int)ceil(tcodeInit + i * codePhaseStep);
        lateCode[i] = caCode[tcode];
      }

      double* promptCode = (double*)malloc(sizeof(double) * blksize);
      tcodeInit = (remCodePhase);
      for (size_t i = 0; i < blksize; i++)
      {
        tcode = (int)ceil(tcodeInit + i * codePhaseStep);
        promptCode[i] = caCode[tcode];

        if (i == (blksize - 1))
        {
          remCodePhase = (tcodeInit + (blksize - 1) * codePhaseStep + codePhaseStep) - 1023.0;
        }
      }

      //// Generate the carrier frequency to mix the signal to baseband---------- -
      double* trigarg = (double*)malloc(sizeof(double) * (blksize + 1));
      for (size_t i = 0; i <= blksize; i++)
      {
        double time = i / settings->samplingFreq;
        trigarg[i] = ((carrFreq * 2.0 * Pi) * time) + remCarrPhase;
      }

      // Get the argument to sin / cos functions
      remCarrPhase = fmod(trigarg[blksize], (2 * Pi));

      // Finally compute the signal to mix the collected data to bandband
      double* iBasebandSignal = (double*)malloc(sizeof(double) * blksize);
      double* qBasebandSignal = (double*)malloc(sizeof(double) * blksize);

      //// Generate the six standard accumulated values-------------------------- -
      // First mix to baseband
      double I_E = 0.0;
      double Q_E = 0.0;
      double I_P = 0.0;
      double Q_P = 0.0;
      double I_L = 0.0;
      double Q_L = 0.0;
#if 0
      FILE* fqBasebandSignal = fopen("qBasebandSignal.txt", "w+");
      FILE* fiBasebandSignal = fopen("iBasebandSignal.txt", "w+");
      fprintf(fqBasebandSignal, "%f\n", lateCode[i]);
      fprintf(fiBasebandSignal, "%f\n", earlyCode[i]);
      fflush(fqBasebandSignal);
      fflush(fiBasebandSignal);
      fclose(fqBasebandSignal);
      fclose(fiBasebandSignal);
#endif
      for (size_t i = 0; i < blksize; i++)
      {
        double qBasebandSignal = cos(trigarg[i]) * (double)rawSignal[i];
        double iBasebandSignal = sin(trigarg[i]) * (double)rawSignal[i];

        I_E += (earlyCode[i] * iBasebandSignal);
        Q_E += (earlyCode[i] * qBasebandSignal);
        I_P += (promptCode[i] * iBasebandSignal);
        Q_P += (promptCode[i] * qBasebandSignal);
        I_L += (lateCode[i] * iBasebandSignal);
        Q_L += (lateCode[i] * qBasebandSignal);
      }

      //// Find PLL error and update carrier NCO----------------------------------

        // Implement carrier loop discriminator(phase detector)
      double carrError = atan(Q_P / I_P) / (2.0 * Pi);

      // Implement carrier loop filter and generate NCO command
      double carrNco = oldCarrNco +
        (tau2carr / tau1carr) * (carrError - oldCarrError) +
        carrError * (PDIcarr / tau1carr);
      oldCarrNco = carrNco;
      oldCarrError = carrError;

      // Modify carrier freq based on NCO command
      carrFreq = carrFreqBasis + carrNco;

      trackResults[channelNr].carrFreq[loopCnt] = carrFreq;

      //// Find DLL error and update code NCO------------------------------------ -
      double codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) /
        (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

      // Implement code loop filter and generate NCO command
      double codeNco = oldCodeNco + (tau2code / tau1code) *
        (codeError - oldCodeError) + codeError * (PDIcode / tau1code);
      oldCodeNco = codeNco;
      oldCodeError = codeError;

      // Modify code freq based on NCO command
      codeFreq = settings->codeFreqBasis - codeNco;

      trackResults[channelNr].codeFreq[loopCnt] = codeFreq;

      //// Record various measures to show in postprocessing----------------------
        // Record sample number(based on 8bit samples)
      trackResults[channelNr].absoluteSample[loopCnt] = fileReadLoc-1;

      trackResults[channelNr].dllDiscr[loopCnt] = codeError;
      trackResults[channelNr].dllDiscrFilt[loopCnt] = codeNco;
      trackResults[channelNr].pllDiscr[loopCnt] = carrError;
      trackResults[channelNr].pllDiscrFilt[loopCnt] = carrNco;

      trackResults[channelNr].I_E[loopCnt] = I_E;
      trackResults[channelNr].I_P[loopCnt] = I_P;
      trackResults[channelNr].I_L[loopCnt] = I_L;
      trackResults[channelNr].Q_E[loopCnt] = Q_E;
      trackResults[channelNr].Q_P[loopCnt] = Q_P;
      trackResults[channelNr].Q_L[loopCnt] = Q_L;

      free(rawSignal);
      free(earlyCode);
      free(lateCode);
      free(promptCode);
      free(trigarg);
      free(iBasebandSignal);
      free(qBasebandSignal);
    }

    // If we got so far, this means that the tracking was successful
    // Now we only copy status, but it can be updated by a lock detector
    // if implemented
    trackResults[channelNr].status = InChannels[channelNr].status;

  } // if a PRN is assigned
}


 /***********************************
 * Public Function Definitions
 ***********************************/

void tracking(
  FILE* fid,
  const Channel_t* InChannels,
  const Settings_t* settings,
  TrackResults_t* trackResults,
  Channel_t* outChannels
)
{

#define FUNC_PREFIX             "tracking(): "
  /*
   Performs code and carrier tracking for all channels.
   Inputs:
   fid - file identifier of the signal record.
    channel - PRN, carrier frequencies and code phases of all
    satellites to be tracked(prepared by preRun.m from
    acquisition results).
    settings - receiver settings.
    Outputs:
   trackResults - tracking results(structure array).Contains
     in - phase prompt outputs and absolute spreading
     code's starting positions, together with other
     observation data from the tracking loops.All are
     saved every millisecond.
    */

  int channelNr = 0;
  // Allocate memory for threads
  thread parallelTrackingThreads[10];

  //// Start processing channels ==============================================
  for (channelNr = 0; channelNr < settings->numberOfChannels; channelNr++)
  {

    /** Call tracking core function in multiple threads
    *
    *
    */
    parallelTrackingThreads[channelNr] = thread(trackingCore,
      channelNr,
      settings,
      fid,
      InChannels,
      trackResults,
      outChannels
    );
  }

  for (channelNr = 0; channelNr < settings->numberOfChannels; channelNr++)
  {
    parallelTrackingThreads[channelNr].join();
  }

}
