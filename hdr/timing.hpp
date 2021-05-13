/**
 * @file timing.hpp
 *
 * @brief C++ header file for defining program execution time and other timing parameters
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

#include <chrono>

#define getTime()                          chrono::steady_clock::now()
#define getElapsedTime(startTime, endTime) chrono::duration_cast<chrono::microseconds>(endTime - startTime).count()


#define Tic1 auto start1 = getTime();
#define Tic2 auto start2 = getTime()
#define Tic3 auto start3 = getTime()
#define Tic4 auto start4 = getTime()
#define Tic5 auto start5 = getTime()
#define Tic6 auto start6 = getTime()

#define Toc1 auto end1 = getTime(); 
#define Toc2 auto end2 = getTime(); 
#define Toc3 auto end3 = getTime(); 
#define Toc4 auto end4 = getTime(); 
#define Toc5 auto end5 = getTime(); 
#define Toc6 auto end6 = getTime(); 

#define GetTime1 getElapsedTime(start1, end1);
#define GetTime2 getElapsedTime(start2, end2);
#define GetTime3 getElapsedTime(start3, end3);
#define GetTime4 getElapsedTime(start4, end4);
#define GetTime5 getElapsedTime(start5, end5);
#define GetTime6 getElapsedTime(start6, end6);

#define Tic(num)  Tic##num
#define Toc(num)  Toc##num
#define Delta(num) GetTime##num;
