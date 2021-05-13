/**
 * @file matcreat.hpp
 *
 * @brief C++ header file for creating SDR MATLAB files
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
#include "mat.h"
#include <stdlib.h>

#define BUFSIZE 256

int writeDoubleValueToMatFile(
  MATFile* pmat,
  double* arr, 
  int row, 
  int col,
  bool isComplex,
  const char* arrName
);

int writeCharValueToMatFile(
  MATFile* pmat,
  const char* arr,
  int len,
  const char* arrName
);


int readStructureFieldArrayDouble(
    MATFile* mfPtr,
    const char* strct,
    const char* field,
    double* outArr,
    int outDim
);

int readStructureFieldArrayINT8(
  MATFile* mfPtr,
  const char* strct,
  const char* field,
  int8_T* outArr,
  int outDim
);
