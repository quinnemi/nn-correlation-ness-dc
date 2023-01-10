#include <boost/multiprecision/cpp_bin_float.hpp>
using namespace boost::multiprecision;
using namespace std;
//enum digit_base_type
//{
//   digit_base_2 = 2,
//   digit_base_10 = 10
//};
//template <unsigned Digits, digit_base_type base = digit_base_10, class Allocator = void, class Exponent = int, ExponentMin = 0, ExponentMax = 0>
//class cpp_bin_float;
//typedef number<cpp_bin_float<128> > float128;

int main() {
    cpp_bin_float_100 b = 2;
    cout << setprecision(numeric_limits<cpp_bin_float_100>::max_digits10)
    << b << endl << exp(b) << endl;
    return 0;
}
