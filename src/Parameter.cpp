
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <concepts> // let's us declare concepts for template constraints

using std::cerr;
using std::endl;
using std::vector;
using std::map;
using std::is_integral_v;
using std::is_floating_point_v;

// Design goals `Parameter`s:
//  - yields samples
//  - has no state (i.e. any state managed by the ABCSMC object)
//  - has no knowledge of the ABCSMC object
//  - does not need to know about other parameters
//
//  challenges to accomplishing this:
//  - need to transform parameters, sometimes in terms of each other
//  - sampling "posterior" or "pseudo" parameters requires state
//
//  way forward:
//  - have "sampling" a posterior / pseudo parameter increment an external state generator (analogous to the RNG for priors)?
//  - have transformations managed by the ABCSMC object?

namespace ABC {

}
