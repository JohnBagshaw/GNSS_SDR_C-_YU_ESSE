/**
 * @file gen_code.cpp
 *
 * @brief code generation source file for SDR functions in C++
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
#include "gen_code.hpp"
#include <math.h>
#include <string.h>
#include <stdlib.h>
 /***********************************
 * Defines
 ***********************************/

 /***********************************
 * Macros
 ***********************************/

 /***********************************
 * Static Variables
 ***********************************/
static const short s_g2Shifts[32] = { 5, 6, 7, 8, 17, 18, 139, 140, 141, 251, 252,
254, 255, 256, 257, 258, 469, 470, 471, 472, 473, 474, 509, 512, 513, 514,
515, 516, 859, 860, 861, 862};


 /***********************************
 * Static Function Definitions
 ***********************************/

static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    }
    else if (u > -0.5) {
      y = u * 0.0;
    }
    else {
      y = ceil(u - 0.5);
    }
  }
  else {
    y = u;
  }

  return y;
}

/***********************************
* Public Function Definitions
***********************************/

void generateCAcode(const Settings_t* settings, double* caCode, int PRN)
{
  int reg1[10] = { 1,1,1,1,1,1,1,1,1,1 };
  int reg2[10] = { 1,1,1,1,1,1,1,1,1,1 };
  int nx;
  int codeValueIndex;
  int saveBit;
  double samplesPerCode;

  int *g1  = new int [1023];
  int *g2  = new int [1023];
  int* g2b = new int [1023];

  /* --- Find number of samples per spreading code ---------------------------- */

  samplesPerCode = settings->samplingFreq / (settings->codeFreqBasis / settings->codeLength);

  if (floor(samplesPerCode) - samplesPerCode)
  {
    samplesPerCode = rt_roundd_snf(samplesPerCode);
  }

  /* --- Generate G1 code ----------------------------------------------------- */
  /* --- Generate all G1 signal chips based on the G1 feedback polynomial ----- */
  for (nx = 0; nx < 1023; nx++) {
    g1[nx] = reg1[9];
    saveBit = reg1[2]^reg1[9];
    for (size_t i = 9; i >= 1; i--)
    {
      reg1[i] = reg1[i-1];
    }
    reg1[0] = saveBit;
  }

  /* --- Generate G2 code ----------------------------------------------------- */
  /* --- Initialize g2 output to speed up the function --- */


  /* --- Generate all G2 signal chips based on the G2 feedback polynomial ----- */
  for (nx = 0; nx < 1023; nx++) {
    g2[nx] = reg2[9];
    saveBit = reg2[1]^reg2[2];
    saveBit ^= reg2[5];
    saveBit ^= reg2[7];
    saveBit ^= reg2[8];
    saveBit ^= reg2[9];
    for (size_t i = 9; i >= 1; i--)
    {
      reg2[i] = reg2[i - 1];
    }
    reg2[0] = saveBit;
  }

  /* --- Shift G2 code -------------------------------------------------------- */
  /* The idea: g2 = concatenate[ g2_right_part, g2_left_part ]; */
  /* --- Form single sample C/A code by multiplying G1 and G2 ----------------- */
  int g2shift = s_g2Shifts[PRN];
  if (g2shift)
  {
    int idx = 0;
    for (nx = 1023 - g2shift; nx < 1023; nx++)
    {
      g2b[idx] = g2[nx];
      idx++;
    }
    for (nx = 0; nx < 1023-g2shift; nx++) 
    {
      g2b[idx] = g2[nx];
      idx++;
    }
  }

  for (nx = 0; nx < 1023; nx++) {
    caCode[nx] = g1[nx] == g2b[nx] ? (double)-1.0 : (double)1.0;
  }

#if 0
  int octalNumber[10] = { 0 };
  int sum = 0;
  int group = 0;
  int digit = 0;
  int idx = 9;

  while (group < 3)
  {
    digit = 0;
    sum = 0;
    while (digit < 3)
    {
      sum += pow(2, digit) * (caCode[idx] == 1 ? 1 : 0);
      idx--;
      digit++;
    }
    octalNumber[3-group] = sum;
    group++;
  }
  octalNumber[3-group] = (caCode[idx] == 1 ? 1 : 0);

  delete[] g1;
  delete[] g2;
  delete[] g2b;
#endif

}

/*
 * Function generates CA codes for all 32 satellites based on the settings
 * provided in the structure "settings". The codes are digitized at the
 * sampling frequency specified in the settings structure.
 * One row in the "caCodesTable" is one C/A code. The row number is the PRN
 * number of the C/A code.
 *
 * caCodesTable = makeCaTable(settings)
 *
 *    Inputs:
 *        settings        - receiver settings
 *    Outputs:
 *        caCodesTable    - an array of arrays (matrix) containing C/A codes
 *                        for all satellite PRN-s
 * Arguments    : const settings *settings
 *                double **caCodesTable
 * Return Type  : void
 */
void makeCaTable(
  const Settings_t* settings, 
  double**caCodesTable,
  int rows, 
  int cols,
  double tS,
  double tC
)
{
  int nx;
  int PRN;
  int codeValueIndex;

  /* --- Prepare the output matrix to speed up function ----------------------- */
  
  /*  C/A chip period in sec */
  /* === For all satellite PRN-s ... */
  codeValueIndex = 0;
  double caCode[1023];
  for (PRN = 0; PRN < rows; PRN++)
  {
    generateCAcode(settings, caCode, PRN);

    /* --- Make index array to read C/A code values ------------------------- */
    /*  The length of the index array depends on the sampling frequency - */
    /*  number of samples per millisecond (because one C/A code period is one */
    /*  millisecond). */

    int count = 1;
    for (nx = 0; nx < cols; nx++) {
      codeValueIndex = ceil((tS * count) / tC);
      if (codeValueIndex > 1023)
      {
        codeValueIndex = 1023;
      }
      caCodesTable[PRN][nx] = caCode[codeValueIndex - 1];
      count++;
    }
    /*  for PRN = 1:32 */
  }
}