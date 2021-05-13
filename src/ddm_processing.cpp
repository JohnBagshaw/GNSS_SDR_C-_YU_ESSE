
#include "ddm_processing.h"

using namespace std;

void DdmProcessing(
  const Settings_t *settings
)
{
#define FUNC_PREFIX "DdmProcessing(): "
    disp("Starting processing...");

    // Define local data structures here
    TrackResults_t trackResults;
    AcqResults_t acqResults;
    WeakAcqSettings_t weakAcqSettings;

    Channel_t inChannels[NUM_CHANNELS];

    FILE* fid = fopen(settings->fileName, "rb");
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
        int16_t* data = new int16_t [BUFFLENGTH]; // 420112
        double* dataDbl = new double[BUFFLENGTH]; // 420112

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

        disp("Doing the demodulation...");
        modulation(dataDbl, settings);

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

          disp("Acquiring satellites...");
          acqResultValid = acquisition(
            dataDbl,
            settings,
            &acqResults
          );
          break;
        }

        case SignalAcquisitionMethod_t::WEAK_SIGNAL_ACQUISITION:
        {

          /** Call normal acquisition function here
          *
          *
          */
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

          break;
        }

        default:
          cout << "No valid choice made for acquisition method" << endl;
          break;

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


        /** Save Acqusition and results to a MAT file for comparison to be done in MATLAB
        *
        *
        */

        // Auto save the acquisitionresults to a file to allow
        // running the positioning solution afterwards.
        disp("Saving DDM Acq results to mat file");
        saveResults(
          &acqResults,
          &trackResults,
          settings,
          true
        );
      }
    }

    // Close the file as processing has finished...
    fclose(fid);
}

static std::mutex mtx;
static std::mutex mtx3;

static void fft(
  fftw_plan p
)
{
  fftw_execute(p);
}

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

static void acquisitionCore(
  const Settings_t* settings,
  AcqResults_t* acqResults,
  int    PRNIdx,
  double* signal1,
  double* signal2,
  double* signal0DC,
  double** caCodeTables,
  bool* sigFound
)
{
#define FUNC_PREFIX " Normal Acquistion (): "
  mtx.lock();
  disp("PRN: %i\n", PRNIdx + 1);
  mtx.unlock();

  // -- - Initialize arrays to speed up the code------------------------------ -
  // Search results of all frequency bins and code shifts(for one satellite)

  // Number of the frequency bins for the given acquisition band(500Hz steps)
  double samplesPerCode = round(settings->samplingFreq / (settings->codeFreqBasis / settings->codeLength));
  int numberOfFrqBins = round(settings->acqSearchBand * 2) + 1;
  double tS = 1 / settings->samplingFreq;
  double tC = 1.0 / settings->codeFreqBasis;
  int fftNumPts = (int)(8 * pow(2, ceil(log(10 * samplesPerCode) / log(2))));
  int uniqFftPts = (int)ceil((fftNumPts + 1) / 2.0);

  double* xCarrier = new double[10 * samplesPerCode];
  double* fftFreqBins = new double[uniqFftPts];
  double** results = new double* [numberOfFrqBins];

  fftw_complex* caCodeFreqDomIn = new fftw_complex[samplesPerCode];
  fftw_complex* caCodeFreqDom = new fftw_complex[samplesPerCode];
  fftw_complex* acqRes1 = new fftw_complex[samplesPerCode];
  fftw_complex* acqRes2 = new fftw_complex[samplesPerCode];
  fftw_complex* xCarrierIn = new fftw_complex[fftNumPts];
  fftw_complex* xCarrierOut = new fftw_complex[fftNumPts];
  fftw_complex** IQfreqDom1In = new fftw_complex * [numberOfFrqBins];
  fftw_complex** IQfreqDom2In = new fftw_complex * [numberOfFrqBins];
  fftw_complex** IQfreqDom1 = new fftw_complex * [numberOfFrqBins];
  fftw_complex** IQfreqDom2 = new fftw_complex * [numberOfFrqBins];
  fftw_complex* IQfreqDom1Loc = new fftw_complex[samplesPerCode];
  fftw_complex* IQfreqDom2Loc = new fftw_complex[samplesPerCode];


  fftw_plan p1 = fftw_plan_dft_1d((int)samplesPerCode, caCodeFreqDomIn, caCodeFreqDom, false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d((int)samplesPerCode, IQfreqDom1Loc, acqRes1, false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p3 = fftw_plan_dft_1d((int)samplesPerCode, IQfreqDom2Loc, acqRes2, false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p4 = fftw_plan_dft_1d(fftNumPts, xCarrierIn, xCarrierOut, true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

  double* frqBinsArr = new double[numberOfFrqBins];

  double scale = 1.0 / samplesPerCode;
  for (int i = 0; i < numberOfFrqBins; ++i)
  {
    results[i] = new double[samplesPerCode];
    IQfreqDom1In[i] = new fftw_complex[samplesPerCode];
    IQfreqDom2In[i] = new fftw_complex[samplesPerCode];
    IQfreqDom1[i] = new fftw_complex[samplesPerCode];
    IQfreqDom2[i] = new fftw_complex[samplesPerCode];

    fftw_plan iqFreqDom1Plans = fftw_plan_dft_1d((int)samplesPerCode, &IQfreqDom1In[i][0], &IQfreqDom1[i][0], true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan iqFreqDom2Plans = fftw_plan_dft_1d((int)samplesPerCode, &IQfreqDom2In[i][0], &IQfreqDom2[i][0], true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);


    double frqBins = settings->iF - (settings->acqSearchBand / 2) * 1000 + 0.5e3 * (i);
    frqBinsArr[i] = frqBins;
    for (size_t j = 0; j < samplesPerCode; j++)
    {
      double sinCar = sin(frqBins * j * 2 * Pi * tS);
      double cosCar = cos(frqBins * j * 2 * Pi * tS);

      // Generate carrier wave frequency grid(0.5kHz step)-----------
      // -- - Generate local sine and cosine------------------------------
      // -- - Convert the baseband signal to frequency domain-------------

      IQfreqDom1In[i][j][0] = sinCar * signal1[j];
      IQfreqDom1In[i][j][1] = cosCar * signal1[j];
      IQfreqDom2In[i][j][0] = sinCar * signal2[j];
      IQfreqDom2In[i][j][1] = cosCar * signal2[j];
    }

    fft(iqFreqDom1Plans);
    fft(iqFreqDom2Plans);

    fftw_destroy_plan(iqFreqDom1Plans);
    fftw_destroy_plan(iqFreqDom2Plans);
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
  fft(p1);

  // Make the correlation for whole frequency band(for all freq.bins)
  for (size_t frqBinIndex = 0; frqBinIndex < numberOfFrqBins; frqBinIndex++)
  {
    // Multiplication in the frequency domain(correlation in time domain)

    for (size_t i = 0; i < samplesPerCode; i++)
    {
      IQfreqDom1Loc[i][0] = (IQfreqDom1[frqBinIndex][i][0] * caCodeFreqDom[i][0]) - (IQfreqDom1[frqBinIndex][i][1] * caCodeFreqDom[i][1]);
      IQfreqDom1Loc[i][1] = (IQfreqDom1[frqBinIndex][i][0] * caCodeFreqDom[i][1]) + (IQfreqDom1[frqBinIndex][i][1] * caCodeFreqDom[i][0]);

      IQfreqDom2Loc[i][0] = (IQfreqDom2[frqBinIndex][i][0] * caCodeFreqDom[i][0]) - (IQfreqDom2[frqBinIndex][i][1] * caCodeFreqDom[i][1]);
      IQfreqDom2Loc[i][1] = (IQfreqDom2[frqBinIndex][i][0] * caCodeFreqDom[i][1]) + (IQfreqDom2[frqBinIndex][i][1] * caCodeFreqDom[i][0]);
    }
    fft(p2);
    fft(p3);

    // -- - Perform inverse DFT and store correlation results------------
    double maxAcqRes1 = 0.0;
    double maxAcqRes2 = 0.0;

    for (size_t i = 0; i < samplesPerCode; i++)
    {
      acqRes1[i][0] *= scale;
      acqRes2[i][0] *= scale;
      acqRes1[i][1] *= scale;
      acqRes2[i][1] *= scale;

      acqRes1[i][0] = acqRes1[i][0] * acqRes1[i][0] + acqRes1[i][1] * acqRes1[i][1];
      acqRes2[i][0] = acqRes2[i][0] * acqRes2[i][0] + acqRes2[i][1] * acqRes2[i][1];

      if (acqRes1[i][0] > maxAcqRes1)
      {
        maxAcqRes1 = acqRes1[i][0];
      }

      if (acqRes2[i][0] > maxAcqRes2)
      {
        maxAcqRes2 = acqRes2[i][0];
      }
    }

    for (size_t i = 0; i < samplesPerCode; i++)
    {
      if (maxAcqRes1 > maxAcqRes2)
      {
        results[frqBinIndex][i] = acqRes1[i][0];
      }
      else
      {
        results[frqBinIndex][i] = acqRes2[i][0];
      }
    }
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

  // Save results
  MATFile* pmat;
  //settings->fileName
  printf("Creating file %s...\n\n", settings->ddmWriteFileName);
  pmat = matOpen(settings->ddmWriteFileName, "w");

  if (pmat == NULL) {
    printf("Error creating file %s\n", settings->ddmWriteFileName);
    printf("(Do you have write permission in this directory?)\n");
  }

  // save acquisition results in the mat file
  int len = numberOfFrqBins * samplesPerCode;
  double* resultsTemp = new double[len];
  for (size_t i = 0; i < numberOfFrqBins; i++)
  {
    for (size_t j = 0; j < samplesPerCode; j++)
    {
      resultsTemp[i * (size_t)samplesPerCode + j] = results[i][j];
    }
  }
  double freqBinIdx = frequencyBinIndex;
  double numFrqBins = numberOfFrqBins;
  double codePhas = codePhase;
  double samplingFreq = settings->samplingFreq;

  writeDoubleValueToMatFile(pmat, resultsTemp, len, 1, false,   "results");
  writeDoubleValueToMatFile(pmat, &freqBinIdx, 1, 1, false,     "frequencyBinIndex");
  writeDoubleValueToMatFile(pmat, &numFrqBins, 1, 1, false,     "numberOfFrqBins");
  writeDoubleValueToMatFile(pmat, &samplesPerCode, 1, 1, false, "samplesPerCode");
  writeDoubleValueToMatFile(pmat, &codePhas, 1, 1, false,       "codePhaseLast");
  writeDoubleValueToMatFile(pmat, &samplingFreq, 1, 1, false,   "samplingFreq");
  writeDoubleValueToMatFile(pmat, &frqBinsArr[0], numberOfFrqBins, 1, false,   "frqBins");
  writeDoubleValueToMatFile(pmat, &*(double*)&settings->acqMethod, 1, 1, false, "weakSigAcqEnabled");

  if (matClose(pmat) != 0) {
    disp("Error closing file %s\n", settings->ddmWriteFileName);
  }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(p3);
  fftw_destroy_plan(p4);

  delete[] resultsTemp;
  delete[] caCodeFreqDomIn;
  delete[] caCodeFreqDom;
  delete[] acqRes1;
  delete[] acqRes2;
  delete[] xCarrierOut;
  delete[] xCarrierIn;
  delete[] fftFreqBins;
  delete[] xCarrier;
  delete[] frqBinsArr;

  //=== Acquisition is over ==================================================

  for (size_t i = 0; i < numberOfFrqBins; i++)
  {
    delete[] results[i];
  }
  delete[] results;

  for (size_t i = 0; i < numberOfFrqBins; i++)
  {
    delete[] IQfreqDom1In[i];
    delete[] IQfreqDom2In[i];

    delete[] IQfreqDom1[i];
    delete[] IQfreqDom2[i];
  }
  delete[] IQfreqDom1In;
  delete[] IQfreqDom2In;
  delete[] IQfreqDom1;
  delete[] IQfreqDom2;
  delete[] IQfreqDom1Loc;
  delete[] IQfreqDom2Loc;

}

static void weakAcquisitionCore(
  int                      PRNIdx,
  const fftw_complex* longSignal,
  const Settings_t* settings,
  const WeakAcqSettings_t* weakAcqSetting,
  WeakAcqResults_t* acqResults,
  double** caCodeTables,
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

  fftw_complex* sigRowFftIn = new fftw_complex[Sblock];
  fftw_complex* sigRowFftOut = new fftw_complex[Sblock];
  fftw_complex* caCodeFftIn = new fftw_complex[Sblock];
  fftw_complex* caCodeFftOut = new fftw_complex[Sblock];
  fftw_complex* sigCaCodeFftIn = new fftw_complex[Sblock];
  fftw_complex* sigCaCodeFftOut = new fftw_complex[Sblock];
  fftw_complex* corrMatrixColIn = new fftw_complex[Nfd];
  fftw_complex* corrMatrixColOut = new fftw_complex[Nfd];
  fftw_complex** corrMatrix = new fftw_complex * [Nfd];
  fftw_complex** corrFftMatrix = new fftw_complex * [Nfd];
  fftw_complex** sigMatrix = new fftw_complex * [Nfd];

  double** caCodeMatrix = new double* [Nfd];

  for (int i = 0; i < Nfd; ++i)
  {
    corrMatrix[i] = new fftw_complex[Sblock];
    corrFftMatrix[i] = new fftw_complex[Sblock];
    sigMatrix[i] = new fftw_complex[Sblock];
    caCodeMatrix[i] = new double[Sblock];
  }

  int PRN = settings->acqSatelliteList[PRNIdx] - 1.0;
  int freqBin = 0;
  int block = 0;
  int blockStart = 0;
  bool foundSignal = false;

  fftw_plan p1 = fftw_plan_dft_1d(Sblock, sigRowFftIn, sigRowFftOut, true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d(Sblock, caCodeFftIn, caCodeFftOut, true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p3 = fftw_plan_dft_1d(Sblock, sigCaCodeFftIn, sigCaCodeFftOut, false ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p4 = fftw_plan_dft_1d(Nfd, corrMatrixColIn, corrMatrixColOut, true ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

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


  // Save results
  MATFile* pmat;

  //settings->fileName
  printf("Creating file %s...\n\n", settings->ddmWriteFileName);
  pmat = matOpen(settings->ddmWriteFileName, "w");

  if (pmat == NULL) {
    printf("Error creating file %s\n", settings->ddmWriteFileName);
    printf("(Do you have write permission in this directory?)\n");
  }

  // save acquisition results in the mat file
  double* resultsTemp = new double[Nfd * Sblock];
  for (size_t i = 0; i < Nfd; i++)
  {
    for (size_t j = 0; j < Sblock; j++)
    {
      resultsTemp[i * (size_t)Sblock + j] = corrFftMatrix[i][j][0];
    }
  }

  double* freqBinArray = new double[Nfd];
  for (size_t i = 0; i < Nfd; i++)
  {
    freqBinArray[i] = settings->iF - (settings->acqSearchBand / 2) * 1000 + 0.5e3 * i;
  }

  double freqBinIdx = freqBin;
  double numFrqBins = Nfd;
  double sBlock = Sblock;
  double codePhas = acqResults->codePhase[PRN];
  double samplingFreq = settings->samplingFreq;
  double weakSigAcqEnabled = settings->acqMethod == SignalAcquisitionMethod_t::WEAK_SIGNAL_ACQUISITION ? 1 : 0;

  writeDoubleValueToMatFile(pmat, resultsTemp,                   Nfd* Sblock, 1, false, "results");
  writeDoubleValueToMatFile(pmat, &freqBinIdx,                   1,           1, false, "frequencyBinIndex");
  writeDoubleValueToMatFile(pmat, &numFrqBins,                   1,           1, false, "numberOfFrqBins");
  writeDoubleValueToMatFile(pmat, &sBlock,                       1,           1, false, "samplesPerCode");
  writeDoubleValueToMatFile(pmat, &codePhas,                     1,           1, false, "codePhaseLast");
  writeDoubleValueToMatFile(pmat, &samplingFreq,                 1,           1, false, "samplingFreq");
  writeDoubleValueToMatFile(pmat, &freqBinArray[0],              Nfd,         1, false, "frqBins");
  writeDoubleValueToMatFile(pmat, &weakSigAcqEnabled,            1,           1, false, "weakSigAcqEnabled");

  delete[] freqBinArray;
  if (matClose(pmat) != 0) {
    disp("Error closing file %s\n", settings->ddmWriteFileName);
  }

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

  delete[] resultsTemp;
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


void DelayDopplerMap(
  const Settings_t* settings
)
{
#define FUNC_PREFIX             "DelayDopplerMap(): "

  AcqResults_t acqResults;

  FILE* fid = fopen(settings->fileName, "rb");
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
      int16_t* data = new int16_t[BUFFLENGTH];
      double* dataDbl = new double[BUFFLENGTH];

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
      disp("Doing the demodulation...");
      modulation(dataDbl, settings);

      switch (settings->acqMethod)
      {
      case SignalAcquisitionMethod_t::NORMAL_SIGNAL_ACQUISITION:
      {
        bool sigFound = false;
        bool sigFoundAny = false;

        // Create two 1msec vectors of data to correlate with and one with zero DC
        double* signal0DC = new double[11 * samplesPerCode];
        double* signal1 = new double[samplesPerCode];
        double* signal2 = new double[samplesPerCode];

        memcpy(signal1, &dataDbl[0], sizeof(double) * samplesPerCode);
        memcpy(signal2, &dataDbl[(int)samplesPerCode], sizeof(double) * samplesPerCode);
        subtractMean(signal0DC, dataDbl, 11 * samplesPerCode);

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

        acquisitionCore(
          settings,
          &acqResults,
          30,
          signal1,
          signal2,
          signal0DC,
          caCodeTables,
          &sigFoundAny
        );

        for (size_t i = 0; i < 32; i++)
        {
          delete[] caCodeTables[i];
        }
        delete[] caCodeTables;
        delete[] signal0DC;
        delete[] signal1;
        delete[] signal2;

        break;
      }
      case SignalAcquisitionMethod_t::WEAK_SIGNAL_ACQUISITION:
      {
        WeakAcqResults_t weakSigAcqResults;
        WeakAcqSettings_t weakAcqSettings;
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

        // Find sampling period
        double tC = 1.0 / settings->codeFreqBasis;

        //  Generate all C/A codes and sample them according to the sampling freq.
        double** caCodeTables = new double* [32];
        for (int i = 0; i < 32; ++i)
        {
          caCodeTables[i] = new double[samplesPerCode];
        }
        makeCaTable(settings, caCodeTables, 32, samplesPerCode, tS, tC);


        disp("(");

        weakAcquisitionCore(
              30,
              baseBandData,
              settings,
              &weakAcqSettings,
              &weakSigAcqResults,
              caCodeTables,
              samplesPerCode
            );
        break;
      }
      default:
      {
        cout << "No suitable acquisition method configured for DDM processing" << endl;
        break;
      }
      }
    }
  }
}