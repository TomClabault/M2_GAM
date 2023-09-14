#ifndef OFFREADER_H
#define OFFREADER_H

#include "mesh.h"

class OffReader
{
public:
    static Mesh read_off(const char* filepath);
};


#endif // OFFREADER_H
