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

class CASim {
public:
    CASim(const char *rule);
    ~CASim();
    
    bool setCells(unsigned int width, unsigned int height, uint8_t max,
                  const uint8_t *cells);
    void getCells(uint8_t *cells);
    
    void step(int n);
    
private:
    struct Rule {
        enum Type {
            LifeLike,
            WireWorld
        };
        
        unsigned short type;
        uint8_t B, S; // parameters for Life-like
    };
    
    Rule rule;
    unsigned int generation;
    unsigned int width;
    unsigned int height;
    uint8_t *cell0;
    uint8_t *cell1; 
};

#endif  // CA_SIM_H
