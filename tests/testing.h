
#ifndef TESTING_H
#define TESTING_H

#include <iostream>

// If macro argument is not true, test is failing
#define IS_TRUE(x) { \
    if (!(x)) { \
        std::cout << __PRETTY_FUNCTION__ << " failed on line " << __LINE__ << std::endl;\
    } else { \
        std::cout << __PRETTY_FUNCTION__ << " passed on line " << __LINE__ << std::endl;\
    } \
}

#endif