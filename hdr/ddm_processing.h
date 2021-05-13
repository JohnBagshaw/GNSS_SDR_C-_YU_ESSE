#pragma once
#include "result_read_write.hpp"
#include "weak_acquisition.hpp"
#include "acquisition.hpp"
#include "gen_code.hpp"
#include "matcreat.hpp"
#include "constants.h"
#include "timing.hpp"
#include "fftw3.h"
#include "mat.h"

#include <mutex>
#include <thread>
#include <math.h>
#include <string.h>
#include <iostream>

void DelayDopplerMap(
  const Settings_t* settings
);

void DdmProcessing(
  const Settings_t* settings
);