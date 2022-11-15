#include <iostream>
#include <vector>

class EnergyLandscape {
    public:
        int L;
        double chemPot;
        std::vector<std::vector<std::vector<double>>> energies;

        EnergyLandscape(int L, std::string pattern) {
            this->L = L;
            this->energies = energies(L, std::vector<std::vector<double>>(L, std::vector<double>(L)));
            if (pattern == "checker") {
                createCheckerPattern();
            }
        }

        void createCheckerPattern(double scale = 1) {
            if (this->L % 2 == 0) {
                throw std::invalid_argument("Length must be odd for checker pattern");
            }

            this->energies
        }
};

int main() {
    std::string a = "hallo";
    if(a == "hallo") {
        std::cout << a << "\n";
    }

    return 0;
}