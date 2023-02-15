#include <iostream>
#include "Euler.h"

int main() {
    std::cout << "ImplicitEuler:\n";
    ImplicitEuler();
    std::cout << "\nExplicitEuler:\n";
    ExplicitEuler();
    std::cout << "\nSymmetric:\n";
    Symmetric();
    std::cout << "\nRungeKutta:\n";
    RungeKutta();

	return 0;
}
