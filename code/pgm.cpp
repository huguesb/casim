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

PGM::PGM()
    : m_owner(true), m_width(0), m_height(0), m_white(0), m_data(0) {
}

PGM::PGM(unsigned int width, unsigned int height, uint8_t white)
    : m_owner(true), m_width(width), m_height(height), m_white(white), m_data(0) {
    m_data = (uint8_t*)malloc(sizeof(uint8_t) * height * width);
}

PGM::PGM(unsigned int width, unsigned int height, uint8_t white, uint8_t *data)
    : m_owner(false), m_width(width), m_height(height), m_white(white), m_data(data) {
    
}

PGM::~PGM() {
    if (m_owner) free(m_data);
}

bool PGM::load(const char *file) {
    if (!m_owner) {
        m_owner = true;
        m_data = 0;
    }
    m_width = 0;
    m_height = 0;
    
    FILE *f = fopen(file, "r");
    if (!f) {
        fprintf(stderr, "Unable to open %s for reading.\n", file);
        return false;
    }
    
    uint8_t magic[2];
    size_t n = fread(magic, sizeof(uint8_t), 2, f);
    if (n != 2 || magic[0] != 'P' || magic[1] != '5') {
        fprintf(stderr, "Invalid PGM file (wrong magic number): %s\n", file);
        return false;
    }
    
    fscanf(f, "%u", &m_width);
    fscanf(f, "%u", &m_height);
    
    unsigned int tmp;
    fscanf(f, "%u", &tmp);
    if (tmp == 0 || tmp > 255) {
        fprintf(stderr, "Invalid PGM file (unsupported white value): %s\n", file);
        m_width = 0;
        m_height = 0;
        return false;
    }
    m_white = tmp;
    
    int c = fgetc(f);
    if (c != ' ' && c != '\n' && c != '\r' && c != '\t') {
        fprintf(stderr, "Invalid PGM file (missing header/data separator): %s\n", file);
        m_width = 0;
        m_height = 0;
        return false;
    }
    
    m_data = (uint8_t*)realloc(m_data, m_width * m_height);
    
    n = fread(m_data, sizeof(uint8_t), m_width * m_height, f);
    if (n != (m_width * m_height)) {
        fprintf(stderr, "Expected %lu bytes of data, got %lu.\n",
                (unsigned long)m_width * (unsigned long)m_height, n);
        m_width = 0;
        m_height = 0;
        return false;
    }
    
    return true;
}

bool PGM::save(const char *file) {
    FILE *f = fopen(file, "w");
    if (!f) {
        fprintf(stderr, "Unable to open %s for writing.\n", file);
        return false;
    }
    if (fprintf(f, "P5\n%u %u\n%u\n", m_width, m_height, m_white) <= 0) {
        fprintf(stderr, "Unable to write PGM header.\n");
        return false;
    }
    size_t sz = (size_t)m_width * (size_t)m_height;
    if (fwrite(m_data, sizeof(uint8_t), sz, f) != sz) {
        fprintf(stderr, "Unable to write PGM data.\n");
        return false;
    }
    return true;
}

void PGM::recode(const std::map<uint8_t, uint8_t>& m) {
    std::map<uint8_t, uint8_t>::const_iterator it, end = m.end();
    unsigned long n = (unsigned long)m_width * (unsigned long)m_height;
    uint8_t *p = m_data;
    while (n) {
        uint8_t v = *p;
        it = m.find(v);
        if (it != end)
            *p = it->second;
        ++p;
        --n;
    }
    it = m.find(m_white);
    if (it != end)
        m_white = it->second;
}

unsigned int PGM::width() const {
    return m_width;
}

unsigned int PGM::height() const {
    return m_height;
}

uint8_t PGM::white() const {
    return m_white;
}

uint8_t* PGM::data() const {
    return m_data;
}
