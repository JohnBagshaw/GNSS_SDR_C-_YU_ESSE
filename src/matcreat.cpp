/**
 * @brief source file for creating SDR MATLAB files in C++
 *
 * Project Title: GNSS-R SDR
 * Author: John Bagshaw
 * Co-author: Surabhi Guruprasad
 * Contact: jotshaw@yorku.ca
 * Supervisor: Prof. Sunil Bisnath
 * Project Manager: Junchan Lee
 * Institution: Lassonde School of Engineering, York University, Canada.

 * MAT-file creation program
 *
 * See the MATLAB External Interfaces/API Guide for compiling information.
 *
 * Calling syntax:
 *
 *   matcreat
 *
 * Create a MAT-file which can be loaded into MATLAB.
 *
 * This program demonstrates the use of the following functions:
 *
 *  matClose
 *  matGetVariable
 *  matOpen
 *  matPutVariable
 *  matPutVariableAsGlobal
 *
 * Copyright 1984-2005 The MathWorks, Inc.
 **/

#include <stdio.h>
#include <string.h> /* For strcmp() */
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include "matcreat.hpp"


int writeDoubleValueToMatFile(
  MATFile* pmat,
  double* arr,
  int row,
  int col,
  bool isComplex,
  const char* arrName
)
{
  int len = row * col;
  int status;

  mxArray* pa;
  double* par, * pai;

  pa = mxCreateDoubleMatrix(row, col, isComplex == true ? mxCOMPLEX : mxREAL);

  if (isComplex == false)
  {
    memcpy((void*)(mxGetPr(pa)), (void*)arr, sizeof(double) * len);
  }
  else
  {
    //mwSize i;
    //par = mxGetPr(pa);
    //pai = mxGetPi(pa);
    //for (i = 0; i < len; i++) {
    //  *par++ = *arr++;
    //  *pai++ = *arr++;
    //}
  }
  status = matPutVariable(pmat, arrName, pa);
  if (status != 0) {
    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
    return(EXIT_FAILURE);
  }

  /* clean up */
  mxDestroyArray(pa);

  return(EXIT_SUCCESS);
}


int writeCharValueToMatFile(
  MATFile* pmat,
  const char* arr,
  int len,
  const char* arrName
) 
{
  int status;
  mxArray* pa;

  pa = mxCreateString(arrName);
  if (pa == NULL) {
    printf("%s :  Out of memory on line %d\n", __FILE__, __LINE__);
    printf("Unable to create string mxArray.\n");
    return(EXIT_FAILURE);
  }
  status = matPutVariable(pmat, arrName, pa);
  if (status != 0) {
    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
    return(EXIT_FAILURE);
  }

  /* clean up */
  mxDestroyArray(pa);

  return(EXIT_SUCCESS);


}

/* Find struct array ARR in MAT-file FILE.
 * Pass field name FIELD to read double arry elements in an array*/
int readStructureFieldArrayINT8(
    MATFile* mfPtr,
    const char* strct,
    const char* field,
    int8_T* outArr,
    int outDim
)
{

  mxArray *pStrct = matGetVariable(mfPtr, strct);
  if (pStrct == NULL) {
    printf("mxArray not found: %s\n", pStrct);
    return(1);
  }

  if (mxGetClassID(pStrct) == mxSTRUCT_CLASS) 
  {
    if (mxGetFieldNumber(pStrct, field) == -1) 
    {
      printf("Field not found: %s\n", field);
    }
    else 
    {
      mwSize nElements;       /* number of elements in array */
      mwIndex eIdx;           /* element index */
      const mxArray* fPtr;    /* field pointer */

      nElements = (mwSize)mxGetNumberOfElements(pStrct);
      for (eIdx = 0; eIdx < nElements; eIdx++)
      {
        if (outDim == eIdx)
        {
          fPtr = mxGetField(pStrct, eIdx, field);
          if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxCHAR_CLASS))
          {
            int i = 0;
            char* ptr = mxArrayToString(fPtr);
            int lElements = (mwSize)mxGetNumberOfElements(fPtr);
            while (i < lElements)
            {
                outArr[i] = ptr[i];
                i++;
            }
          }
        }
      }
    }
  }
  else {
    printf("%s is not a structure\n", field);
  }
  mxDestroyArray(pStrct);
  return(0);
}

int readStructureFieldArrayDouble(
    MATFile* mfPtr,
    const char* strct,
    const char* field,
    double* outArr,
    int outDim
)
{

    mxArray* pStrct = matGetVariable(mfPtr, strct);
    if (pStrct == NULL) {
        printf("mxArray not found: %s\n", pStrct);
        return(1);
    }

    if (mxGetClassID(pStrct) == mxSTRUCT_CLASS)
    {
        if (mxGetFieldNumber(pStrct, field) == -1)
        {
            printf("Field not found: %s\n", field);
        }
        else
        {
            mwSize nElements;       /* number of elements in array */
            mwIndex eIdx;           /* element index */
            const mxArray* fPtr;    /* field pointer */

            nElements = (mwSize)mxGetNumberOfElements(pStrct);
			for (eIdx = 0; eIdx < nElements; eIdx++) 
            {
                if (outDim == eIdx)
                {
                    fPtr = mxGetField(pStrct, eIdx, field);
                    if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS))
                    {
                        int i = 0;
                        double* ptr = mxGetPr(fPtr);
                        int lElements = (mwSize)mxGetNumberOfElements(fPtr);
                        while (i < lElements)
                        {
                            outArr[i] = ptr[i];
                            i++;
                        }
                    }
                }
			}
        }
    }
    else {
        printf("%s is not a structure\n", field);
    }
    mxDestroyArray(pStrct);
    return(0);
}
