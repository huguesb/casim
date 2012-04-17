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

#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage : %s <input.pgm> <output.pgm> <oldCol0:newCol1>+\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    
    PGM pgm;
    std::map<uint8_t, uint8_t> colormap;
    
    for (int i = 3; i < argc; ++i) {
        unsigned int oval, nval;
        sscanf(argv[i], "%u:%u", &oval, &nval);
        colormap[oval & 0xff] = nval & 0xff;
    }
    
    if (!pgm.load(argv[1])) return EXIT_FAILURE;
    pgm.recode(colormap);
    if (!pgm.save(argv[2])) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
