#include <iostream>
#include <cassert>

int add(int a, int b) {
    return a + b;
}

int main() {
    int result = add(2, 3);
    assert(result == 5);
    std::cout << "Test passed\n";
    return 0;
}