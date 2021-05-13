/**
 * @file structures.hpp
 *
 * @brief C++ header file for SDR structures
 *
 * Project Title: GNSS-R SDR
 * Author: John Bagshaw
 * Co-author: Surabhi Guruprasad
 * Contact: jotshaw@yorku.ca
 * Supervisor: Prof. Sunil Bisnath
 * Project Manager: Junchan Lee
 * Institution: Lassonde School of Engineering, York University, Canada.
 **/

#pragma once


#include "constants.h"

typedef enum DataType_t {
  DATA_TYPE_INT8 = 0,
  DATA_TYPE_INT16
};

typedef enum class SignalAcquisitionMethod_t {
  NORMAL_SIGNAL_ACQUISITION = 0,
  WEAK_SIGNAL_ACQUISITION = 1,
  HALF_BIT_METHOD    = 2
}SignalAcquisitionMethod_t;

typedef struct {
  double E;
  double N;
  double U;
} TruePosition_t;

typedef struct {
  int gnssStart;
  double msToProcess;
  int numberOfChannels;
  double skipNumberOfBytes;
  char   fileName      [MAX_FILE_NAME_LEN_CHARS];
  char   readFileName  [MAX_FILE_NAME_LEN_CHARS];
  char   writeFileName [MAX_FILE_NAME_LEN_CHARS];
  char   ddmWriteFileName[MAX_FILE_NAME_LEN_CHARS];
  DataType_t dataType;
  double iF;
  double samplingFreq;
  double codeFreqBasis;
  double codeLength;
  double skipAcquisition;
  double acqSatelliteList[32];
  double acqSearchBand;
  double acqThreshold;
  double dllDampingRatio;
  double dllNoiseBandwidth;
  double dllCorrelatorSpacing;
  double pllDampingRatio;
  double pllNoiseBandwidth;
  double navSolPeriod;
  double elevationMask;
  double useTropCorr;
  bool plotTracking;
  TruePosition_t truePosition;
  double c;
  double startOffset;
  SignalAcquisitionMethod_t acqMethod;
  double acqMS;

} Settings_t;

typedef struct WeakAcqSettings_s
{

  double PIT;
  double totalSamples;

  // TODO: add more here ??

}WeakAcqSettings_t;

typedef struct {
  double carrFreq[32];
  double codePhase[32];
  double peakMetric[32];
} AcqResults_t;

typedef struct {
  double carrFreq[32];
  double codePhase[32];
  double peakMetric[32];
} WeakAcqResults_t;

typedef struct {
  bool isMatEqual_carrFreq;
  bool isMatEqual_codePhase;
  bool isMatEqual_peakMetric;
} IsMatEqualAcqResults_t;


typedef struct {
  double PRN;
  double acquiredFreq;
  double codePhase;
  char status;
} Channel_t;

typedef struct {

  char   status;
  double PRN;
  double* absoluteSample;
  double* codeFreq;
  double* carrFreq;
  double* I_P;
  double* I_E;
  double* I_L;
  double* Q_E;
  double* Q_P;
  double* Q_L;
  double* dllDiscr;
  double* dllDiscrFilt;
  double* pllDiscr;
  double* pllDiscrFilt;

} TrackResults_t;

typedef struct {
  bool isMatEqual_status;
  bool isMatEqual_PRN;
  bool isMatEqual_absoluteSample;
  bool isMatEqual_codeFreq;
  bool isMatEqual_carrFreq;
  bool isMatEqual_I_P;
  bool isMatEqual_I_E;
  bool isMatEqual_I_L;
  bool isMatEqual_Q_E;
  bool isMatEqual_Q_P;
  bool isMatEqual_Q_L;
  bool isMatEqual_dllDiscr;
  bool isMatEqual_dllDiscrFilt;
  bool isMatEqual_pllDiscr;
  bool isMatEqual_pllDiscrFilt;
} IsMatEqualTrackResults_t;
