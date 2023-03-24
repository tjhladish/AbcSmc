#include <iostream>
#include <vector>

#include <stdlib.h>
#include <time.h>

using namespace std;

class Par {
  public:
    Par(): _vmin(0), _vmax(0), _state(0), _is_random(false) {}
    Par(int vmin, int vmax, int state, bool is_random): _vmin(vmin), _vmax(vmax), _state(state), _is_random(is_random) {}

    int min() { return _vmin; }
    void set_min(int val) { _vmin = val; }
    int max() { return _vmax; }
    void set_max(int val) { _vmax = val; }
    int state() { return _state; }
    void set_state(int val) { _state = val; }
    void increment() { ++_state; }
    void reset() { _state = _vmin; }
    bool is_random() { return _is_random; }

  private:
    int _vmin, _vmax, _state;
    bool _is_random;
};

void report_state(const vector<Par*> &pars) {
    for (auto p: pars) {
        cout << p->state() << " ";
    }
    cout << endl;
}

void next_combination(vector<Par*> &pars) {
    bool increment_nonrandom_par = true;
    for (auto p: pars) {
        if (p->is_random()) {
            p->set_state(p->min() + (rand() % (p->max() - p->min() + 1)));
        } else {
            if (increment_nonrandom_par) {
                if (p->state() == p->max()) {
                    p->reset();
                } else {
                    p->increment();
                    increment_nonrandom_par = false;
                }
            }
        }
    }
}

int main() {
        // constructor args: min_val, max_val, starting state, whether to sample randomly (true) or increment (false)
    Par* p1 = new Par(0, 9, 0, false);
    Par* p2 = new Par(3, 5, 3, true);
    Par* p3 = new Par(-3, 5, 3, false);
    vector<Par*> pars = {p1, p2, p3};

    for (int i = 0; i < 200; ++i) {
        report_state(pars);
        next_combination(pars);
    }
}
