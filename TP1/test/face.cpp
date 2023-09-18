#include "face.h"

std::ostream& operator << (std::ostream& os, const Face& face)
{
    os << "Face[" << face.m_a << ", " << face.m_b << ", " << face.m_c << "]";

    return os;
}
