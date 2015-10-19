#ifndef MAT_CONFIG_HPP
#define MAT_CONFIG_HPP

/*may need to replace with correct alloca header*/
#include <stdlib.h>
#include <alloca.h>
#include <assert.h>
//stack allocate
#define MAT_ALLOCA(sz) alloca(sz)
//on false, abort. or throw.
#define MAT_ASSERT(expr) assert(expr)
//default scalar type
#define MAT_SCALAR_TYPE float
//#define MAT_IPML_NO_MATH // for platforms without math.h

#endif
