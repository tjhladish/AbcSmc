#ifndef ENUM_MACROS_H
#define ENUM_MACROS_H

// provide pieces for stringify-ing Event<...>
#include <string>
#include <iostream>

// defining these macros to reduce language noise
// TODO: do this the "right" way i.e. port FOR_EACH ala https://www.scs.stanford.edu/~dm/blog/va-opt.html
// so that MAKE_ENUM can also include defining << operator, other conveniences

// templating tools to enable convenient "EnumType e = from_string(str);" syntax
template<class ET>
inline ET hidden_fs(const std::string &/* */) { return new ET; }

struct SEWrapper {
    const std::string str;
    SEWrapper(const std::string str) : str(str) {}
    template<typename T>
    operator T () const { return hidden_fs<T>(str); } 
};

inline SEWrapper from_string(const std::string str) {
    return SEWrapper(str);
};

// now various macros to build up enums conveniently

// fundamental tools
#define CONCAT(a, b) CONCAT_(a, b)
#define CONCAT_(a, b) a ## b
#define COMMA(X) X,
#define QUOTE(X) #X

// these pieces are the typical Enum => String and vice versa bits
#define SSTR(X) case X: return QUOTE(X); break;
#define STRTOE(X) if (estr == QUOTE(X)) { return X; } else
// if using "and counter" enum style - i.e. last enum of an EnumType is N_EnumType
#define N_ENUM(ET) CONCAT(N_,ET)

// constructs an enum; the VA args should COMMA(Enum1) COMMA(Enum2) etc
#define MAKE_ENUM(ET, ...) enum ET { \
  __VA_ARGS__ \
  N_ENUM(ET) \
};

#define STRINGIFY_ENUM(ET, switchblock)\
inline std::string to_string(const ET &e) { \
  switch (e) {\
    switchblock\
    default: exit(-1); break;\
  }\
};\
\
inline std::ostream& operator<<(std::ostream &out, const ET &e) { \
  return out << to_string(e); \
};

#define PARSE_ENUM(ET, ifblock)\
template<>\
inline ET hidden_fs<ET>(const std::string &estr) { \
  ifblock \
  STRTOE(CONCAT(N_,ET)) { \
    std::cerr << "Failed to parse " << #ET << " from: " << estr << " - exiting." << std::endl; \
    exit(-1); \
  }\
};

#define CONSTRUCTENUM(ENUMM)\
ENUMM(MAKE_ENUM, COMMA) \
ENUMM(STRINGIFY_ENUM, SSTR) \
ENUMM(PARSE_ENUM, STRTOE)

/* example usage - think of this as replacing normal Enum declaration
#define ENUMDEMO(MACRO, SUBM) MACRO(EnumDemo,SUBM(FOO) SUBM(BAR))
CONSTRUCTENUM(ENUMDEMO)
*/

#endif