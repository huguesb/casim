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

#ifndef CA_SIM_H
#define CA_SIM_H

#include <stdint.h>
#include <stdlib.h>

struct CARule;

class CASim {
public:
    CASim(const char *rule);
    ~CASim();
    
    bool setCells(unsigned int width, unsigned int height, uint8_t max,
                  const uint8_t *cells);
    void getCells(uint8_t *cells);
    
    void step(int n);
    
private:
    CARule *rule;
    unsigned int generation;
    unsigned int width;
    unsigned int height;
    size_t pitch;
    uint8_t *cell0;
    uint8_t *cell1;
    uint32_t *work0;
    uint32_t *work1;
    uint32_t *workOffset;
};

#endif  // CA_SIM_H
