#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>


template <typename T>
inline std::string toString (const T& t) {
    std::stringstream ss;
    ss  << std::setprecision(30) << t;
    return ss.str();
}


int main() 
{
    //float f = 23.12345678900000000010000000001000000000100000000010000000001;
    float f = (long double) (23.0123456789012345678901234567890);
    std::string f_str = std::to_string(f);
    std::cout << f_str.length() << " " << f_str << '\n';

    //double g = 23.12345678900000000010000000001000000000100000000010000000001;
    double g = (long double) (23.0123456789012345678901234567890L);
    std::string g_str = std::to_string(g);
    std::cout << g_str.length() << " " << g_str << '\n';

    //long double h = 23.12345678900000000010000000001000000000100000000010000000001;
    long double h = (long double) (23.0123456789012345678901234567890L);
    std::string h_str = std::to_string(h);
    std::cout << h_str.length() << " " << h_str << '\n';

    f_str = toString(f);
    std::cout << f_str.length() << " " << f_str << '\n';

    g_str = toString(g);
    std::cout << g_str.length() << " " << g_str << '\n';

    h_str = toString(h);
    std::cout << h_str.length() << " " << h_str << '\n';

    for (int i = 0; i<50; i++) {
        double f = i - 0.00000000000001;
        int j = (int) f;
        std::cout << std::setprecision(20) << i << " " << f << " " << j << std::endl;
        std::cout << std::setprecision(8) << i << " " << f << " " << j << std::endl;
    }


}
