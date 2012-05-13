/****************************************************************************
** Copyright (c) 2012 Hugues Bruant <hugues@cmu.edu>
** All rights reserved.
**
** This file is part of a school project and licensed under the terms of FreeBSD
** license (2-clause BSD also refered to as Simplified BSD License)
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
****************************************************************************/

#include "pgm.h"
#include "casim.h"
#include "cycleTimer.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage : %s <input.pgm> <output.pgm> <rule> <steps>\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    
    double ref = CycleTimer::currentSeconds();
    
    PGM input;
    CASim sim(argv[3]);
    
    if (!input.loadPattern(argv[1])) { return EXIT_FAILURE; }
    
    double ldr = CycleTimer::currentSeconds();
    
    sim.setCells(input.width(), input.height(), input.white(), input.data());
    
    double set = CycleTimer::currentSeconds();
    
    sim.step(atoi(argv[4]));
    
    double stp = CycleTimer::currentSeconds();
    
    sim.getCells(input.data());
    
    double get = CycleTimer::currentSeconds();
    
    if (!input.save(argv[2])) { return EXIT_FAILURE; }
    
    double str = CycleTimer::currentSeconds();
    
    fprintf(stderr, "Load  : %lf ms\n", (ldr - ref)*1000.0);
    fprintf(stderr, "Set   : %lf ms\n", (set - ldr)*1000.0);
    fprintf(stderr, "Sim   : %lf ms\n", (stp - set)*1000.0);
    fprintf(stderr, "Get   : %lf ms\n", (get - stp)*1000.0);
    fprintf(stderr, "Store : %lf ms\n", (str - get)*1000.0);
    fprintf(stderr, "Total : %lf ms\n", (str - ref)*1000.0);
    
    return EXIT_SUCCESS;
}
