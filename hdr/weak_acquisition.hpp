/**
 * @file weak_sig_acquisition.hpp
 *
 * @brief C++ header file for SDR weak signal acquisition
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

 /***********************************
 * Includes
 ***********************************/
#include "fftw3.h"
#include "structures.hpp"
#include <stdio.h>

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

 /***********************************
 * Public Function Declarations
 ***********************************/


bool weak_acquisition(
  const fftw_complex* longSignal,
  const Settings_t* settings,
  const WeakAcqSettings_t* weakAcqSetting,
  WeakAcqResults_t* acqResults
);





