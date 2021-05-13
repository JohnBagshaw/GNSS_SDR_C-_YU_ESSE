/**
 * @file constants.h
 *
 * @brief C++ header file for SDR constants
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

#define _USE_MATH_DEFINES

#define NUM_CHANNELS            (8)
#define MAX_FILE_NAME_LEN_CHARS (200)

// Control how many milliseconds of code to combine in normal acquisition.
// This value must be atleast two.
#define BUFFLENGTH              (1600000)
#define Pi                      ( M_PI )//pi value

