#include <cassert>
#include <iostream>

int main() {
    std::cout << "[smoke] Running..." << std::endl;
    assert(1 + 1 == 2);
    std::cout << "[smoke] OK" << std::endl;
    return 0;
}
