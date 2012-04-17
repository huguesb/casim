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

#ifndef PGM_H
#define PGM_H

#include <stdint.h>

class PGM {
public:
    PGM();
    PGM(unsigned int width, unsigned int height, uint8_t white = 255);
    PGM(unsigned int width, unsigned int height, uint8_t white, uint8_t *data);
    ~PGM();
    
    bool load(const char *file);
    bool save(const char *file);
    
    unsigned int width() const;
    unsigned int height() const;
    uint8_t white() const;
    uint8_t* data() const;
    
private:
    bool m_owner;
    unsigned int m_width;
    unsigned int m_height;
    uint8_t m_white;
    uint8_t *m_data;
};

#endif  // PGM_H
