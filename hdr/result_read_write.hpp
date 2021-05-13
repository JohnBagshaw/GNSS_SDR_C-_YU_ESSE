/**
  * @file result_read_write.hpp
  *
  * @brief C++ header file for reading/writing SDR results
  *
  * Project Title: GNSS-R SDR
  * Author: John Bagshaw
  * Co-author: Surabhi Guruprasad
  * Contact: jotshaw@yorku.ca
  * Supervisors: Prof. Sunil Bisnath
  * Project Manager: Junchan Lee
  * Institution: Lassonde School of Engineering, York University, Canada.
  **/

#pragma once

 /***********************************
 * Includes
 ***********************************/
#include "structures.hpp"
#include <chrono>
#include <stdio.h>

using namespace std;
 /***********************************
 * Defines
 ***********************************/
//#define DONT_DO_PRINTF

 /***********************************
 * Macros
 ***********************************/
#ifndef DONT_DO_PRINTF
#define disp(f_, ...) printf(FUNC_PREFIX); printf((f_), __VA_ARGS__); printf("\n") 
#else 
#define disp(f_, ...)
#endif

 /***********************************
 * Static Variables
 ***********************************/

 /***********************************
 * Static Function Definitions
 ***********************************/

 /***********************************
 * Public Function Declarations
 ***********************************/

void preRun(
  AcqResults_t* acqResults,
  const Settings_t* settings,
  Channel_t* channel
);

void showChannelStatus(
  const Settings_t* settings,
  const Channel_t* channel
);

void saveResults(
  AcqResults_t* acqResults,
  TrackResults_t* trackResults,
  const Settings_t* settings,
  bool disableTrackResultSave
);

void compareAcquisitionResults(
  const Settings_t* settings,
  AcqResults_t* acqResults,
  IsMatEqualAcqResults_t* isMatEqualAcqResults
);

void compareTrackingResults(
  TrackResults_t* trackResults,
  IsMatEqualTrackResults_t* isMatEqualTrackResults,
  const Settings_t* settings
);