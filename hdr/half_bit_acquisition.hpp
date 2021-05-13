/**
 * @file half_bit_acquisition.hpp
 *
 * @brief C++ header file for alternate half-bit acquisition method
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

bool half_bit_acquisition(
  const double* longSignal,
  const Settings_t* settings,
  AcqResults_t* acqResults,
  int dataLength
);





