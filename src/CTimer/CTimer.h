/**
  Copyright (c) 2012-2015 "Bordeaux INP, Bertrand LE GAL"
  [http://legal.vvv.enseirb-matmeca.fr]

  This file is part of LDPC_C_Simulator.

  LDPC_C_Simulator is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CLASS_CTimer
#define CLASS_CTimer

#define USE_CUSTOM_TIMER

#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

#ifndef USE_CUSTOM_TIMER
    #include "boost/timer/timer.hpp"
    using boost::timer::cpu_timer;
    using boost::timer::cpu_times;
    using boost::timer::nanosecond_type;
#else
    #include <time.h>
    #include <sys/time.h>
    #include <stdlib.h>

    #ifdef __MACH__
        #include <mach/clock.h>
        #include <mach/mach.h>
    #endif
#endif


class CTimer
{
    
private:
#ifndef USE_CUSTOM_TIMER
    cpu_timer timer;
#else
    timespec t_start;
    timespec t_stop;
#endif
    bool isRunning;
    
protected:
#ifdef USE_CUSTOM_TIMER
    long diff_ns(timespec start, timespec end);
    long diff_sec(timespec start, timespec end);
#endif

public:
    CTimer(bool _start);
    CTimer();
    void start();
    void stop();
    void reset();
    long get_time_ns();
    long get_time_us();
    long get_time_ms();
    long get_time_sec();
};

#endif
