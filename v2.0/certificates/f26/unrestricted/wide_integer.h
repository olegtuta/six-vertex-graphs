#if defined(_MSC_VER) && !defined(__clang__)
#include <__msvc_int128.hpp>
using WideInteger = std::_Signed128;
#else
using WideInteger = __int128;
#endif
