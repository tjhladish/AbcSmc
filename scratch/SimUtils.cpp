#include <compare>
#include <SimUtils.h>

// Calculates the number of times the data series crosses the median
float_type median_crossings(const Col & data, const float_type m) {
    if (data.size() < 2) return 0.0;
    
    size_t mc = 0;

    // current and next are like cursors that trace the data series
    // this uses the <=> operator: (a <=> b) returns equivalent if a == b, less if a < b, and greater if a > b
    partial_ordering current = data[0] <=> m;
    if (current == partial_ordering::equivalent) mc++; // increment if we're starting at the median
    for (auto pos : data) { // for each data point (n.b. the first pass here is a no-op: next == current)
        partial_ordering next = pos <=> m;
        if ((next != current) and (current != partial_ordering::equivalent)) mc++; // just crossed or touched the median
        current = next;
    }
    return (static_cast<float_type>(mc)/(data.size() - 1.0));
}

float_type median_crossings(const Col & data) {
    if (data.size() < 2) {
        return 0;
    } else {
        return median_crossings(data, median(data));
    }
}