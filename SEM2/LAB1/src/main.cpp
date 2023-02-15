#include <iostream>
#include "Euler.h"

int main() {
    std::cout << "ImplicitEuler: ";
    ImplicitEuler();
    std::cout << "\nExplicitEuler:\n";
    ExplicitEuler();

	return 0;
}
