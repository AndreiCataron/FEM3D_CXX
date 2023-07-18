#include <iostream>

int main() {
    std::vector<int> a = {1, 2, 3};
    std::vector<int> b = std::vector<int>(a.begin() + 2, a.begin() + 3);

    for(auto i : b){
        std::cout << i << ' ';
    }
}