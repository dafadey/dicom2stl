#pragma once

#include <vector>
#include <string>
typedef unsigned char BYTE;

namespace base64 {

std::string encode(BYTE const* buf, unsigned int bufLen);
std::vector<BYTE> decode(std::string const&);
}