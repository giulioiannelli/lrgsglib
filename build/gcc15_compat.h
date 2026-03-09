/*
 * gcc15_compat.h — Workaround for GCC 15 + old conda sysroot
 *
 * GCC 15's <ctime> does `using ::timespec_get;` when
 * _GLIBCXX_HAVE_TIMESPEC_GET is set in bits/c++config.h.
 * Old conda sysroots (glibc < 2.28) don't provide timespec_get
 * in <time.h>, causing a compilation error.
 *
 * Fix: force-include this header via -include to undefine the macro
 * AFTER c++config.h sets it but BEFORE <ctime> uses it.
 *
 * Usage:  CXXFLAGS += -include /path/to/gcc15_compat.h
 */
#ifdef __cplusplus
#include <bits/c++config.h>
#ifdef _GLIBCXX_HAVE_TIMESPEC_GET
#undef _GLIBCXX_HAVE_TIMESPEC_GET
#endif
#endif
