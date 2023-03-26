
#include "testing.h"
#include "pls.h"

// Function to test
bool test_order(const Col &test, const std::vector<size_t> &target) {
    auto res = ordered(test);
    bool same = true;
    for (size_t i = 0; (i < test.size()) and same; i++) {
        same = same and (res[i] == target[i]);
    }
    return same;
}

void series_order() {
    Row ref = Row(3);
    ref << 1, 2, 3;
    std::vector<size_t> ord = { 0, 1, 2 };
    IS_TRUE(test_order(ref, ord));
    ref << 2, 1, 3;
    ord = { 1, 0, 2 };
    IS_TRUE(test_order(ref, ord));
    // TODO more test cases?
}

int main(void) {
    series_order();
}