#include <iostream>
#include <omp.h>

class A {
public:
    virtual void saysth() = 0;
};

class B : public A {
private:
    int x;
};

class C : public B {
public:
    void saysth() {
        std::cout << "Salut";
    }
};

int main() {
    std::vector<int> a = {1, 2, 3};
    std::vector<int> b = std::vector<int>(a.begin() + 2, a.begin() + 3);

    for(auto i : b){
        std::cout << i << ' ';
    }

    C obj;
    obj.saysth();

    #pragma omp parallel for default(none)

    for (int i = 0; i < 10000000; i++) {
        int a = 0;
        a = a * a - a;
        std::vector<int> v;
        v.push_back(a);
    }
}