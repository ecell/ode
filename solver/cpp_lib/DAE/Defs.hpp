#ifndef __DEFS_HPP
#define __DEFS_HPP

#include <boost/call_traits.hpp>
#include <boost/smart_ptr.hpp>

#include <vector>
#include <limits>

/**
    Declare a new type and its pointer,
    const pointer, reference, and const reference types. For example
    DECLARE_TYPE(Dword, VeryLongTime);
    @param mydecl The base type
    @param mytype The new type
*/

#define DECLARE_TYPE(mydecl, mytype)  \
typedef mydecl         mytype;         \
typedef mytype *       mytype ## Ptr;  \
typedef const mytype * mytype ## Cptr; \
typedef mytype &       mytype ## Ref;  \
typedef const mytype & mytype ## Cref;

// Types

template <typename T>
class Param
{
public:
    typedef typename boost::call_traits<T>::param_type type;
};

typedef double Real;
//typedef Param<Real>::type RealParam;
typedef std::vector<Real> VariableArray;

typedef Real Time;
typedef Param<Time>::type TimeParam;

typedef long int Integer;
typedef unsigned long int UnsignedInteger;

typedef std::vector<Real> RealVector;

//20100224追加
struct StatusEvent
{
    Integer variable_index;
    bool variable_flag;//0ならDifferential, 1ならAlgebraic
    Real threshold;
    Integer status_code;
};
typedef std::vector<StatusEvent> StatusEventArray;

//! Infinity.
const Real INF(std::numeric_limits<double>::infinity());

// *******************************************
// Define the NULLPTR
// *******************************************

const int NULLPTR(0);



#endif /* __DEFS_HPP */
