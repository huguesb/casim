/****************************************************************************
** Copyright (c) 2012 Hugues Bruant <hugues@cmu.edu>
** All rights reserved.
**
** This file may be used under the terms of the GNU General Public License
** version 3 as published by the Free Software Foundation.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
****************************************************************************/

#include "pgm.h"
#include "casim.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage : %s <input.pgm> <output.pgm> <rule> <steps>\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    
    PGM input;
    CASim sim(argv[3]);
    
    if (!input.load(argv[1])) { return EXIT_FAILURE; }
    sim.setCells(input.width(), input.height(), input.white(), input.data());
    
    sim.step(atoi(argv[4]));
    
    sim.getCells(input.data());
    if (!input.save(argv[2])) { return EXIT_FAILURE; }
    return EXIT_SUCCESS;
}
