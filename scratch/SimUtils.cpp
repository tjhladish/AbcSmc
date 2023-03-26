#include <compare>
#include <SimUtils.h>
#include <AbcSmc/AbcUtil.h> // for median, Col

int _mc_pos(const float_type v, const float_type m) {
    enum Position {ABOVE, BELOW, AT};
    Position p;
    if (v > m) { p = ABOVE; }
    else if (v < m) { p = BELOW; }
    else { p = AT; }
    return (int) p;
}

// Calculates the number of times the data series crosses the median
float_type median_crossings(const Col &data, const float_type m) {
    int mc = 0;
    if (data.size() < 2) return mc;

    enum Position {ABOVE, BELOW, AT};
    // current and next are like cursors that trace the data series
    Position current, next;
    current = (Position) _mc_pos(*(data.begin()), m);
    if (current == AT) mc++; // increment if we're starting at the median
    for (T* it = (data.begin()+1); it != data.end(); it++) {
        next = (Position) _mc_pos(*it, m);
        if (next != current and current != AT) mc++; // just crossed or touched the median
        current = next;
    }
    return ((float_type) mc)/(data.size()-1);
}

float_type median_crossings(const Col &data) {
    if (data.size() < 2) {
        return 0;
    } else {
        return median_crossings(data, median(data));
    }
}


// Calculates the number of times the data series crosses the median
// float_type median_crossings(const Col &data, const float_type m) {
//     if (data.size() < 2) return 0.0;
    
//     size_t mc = 0;

//     // current and next are like cursors that trace the data series
//     // this uses the <=> operator: (a <=> b) returns equivalent if a == b, less if a < b, and greater if a > b
//     partial_ordering current = data[0] <=> m;
//     if (current == partial_ordering::equivalent) mc++; // increment if we're starting at the median
//     for (auto pos : data) { // for each data point (n.b. the first pass here is a no-op: next == current)
//         partial_ordering next = pos <=> m;
//         if ((next != current) and (current != partial_ordering::equivalent)) mc++; // just crossed or touched the median
//         current = next;
//     }
//     return (static_cast<float_type>(mc)/(data.size() - 1.0));
// }

// float_type median_crossings(const Col &data) {
//     if (data.size() < 2) {
//         return 0;
//     } else {
//         return median_crossings(data, median(data));
//     }
// }