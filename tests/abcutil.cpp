
#include "testing.h"
#include "AbcUtil.h"

// Function to test
bool test_colwise_z_scores(const Mat2D &test, const Mat2D &target) {
    Mat2D res = PLS::colwise_z_scores(test);
    return (target - res).array().square().sum() < 1e-6;
}

void series_colwise_z_scores() {
    Mat2D ref = Mat2D(3, 3), stand = Mat2D(3, 3);
    ref << 1, 1, 1,
           2, 3, 4,
           3, 5, 7;
    stand << -1, -1, -1,
             0, 0, 0,
             1, 1, 1;
    IS_TRUE(test_colwise_z_scores(ref, stand));
    // TODO more test cases?
}

bool test_euclidean(const Mat2D &test, const Row &target, const Col &distref) {
    Col res = ABC::euclidean(test, target);
    return (res - distref).norm() < 1e-6;
}

void series_euclidean() {
    Mat2D ref = Mat2D(2, 2);
    Row tar = Row(2);
    Col distref = Col(2);
    ref << 1, 1,
           3, 3;
    tar << 1, 1;
    distref << 0,
               2.828427;
               
    IS_TRUE(test_euclidean(ref, tar, distref));
    // TODO more test cases?
}

int main(void) {
    series_colwise_z_scores();
    series_euclidean();
}