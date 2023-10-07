#include <iostream>
//#include <omp.h>

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

struct s {
    int x;

    s(int a) {
        x = a;
    }
};

struct ss : s {
    int y;
    ss(int a, int b) : s(a){
        y = b;
    }
};

struct parent {
    int x = 0;
    explicit parent(int a) {
        x = a;
    }
};

struct other {
    int y;
    other(int a) {
        y = a;
    }
};

class D {
public:
    void printsth() {std::cout << "SALUT";}
};

class E : public D {
private:
    int x;
public:
    explicit E(int a){x=a;}
};

int main() {
//    std::vector<int> a = {1, 2, 3};
//    std::vector<int> b = std::vector<int>(a.begin() + 2, a.begin() + 3);
//
//    for(auto i : b){
//        std::cout << i << ' ';
//    }
//
//    C obj;
//    obj.saysth();
//
//    s model = {5};
//
//    ss model2 = {6, 7};
//    std::cout << '\n' << model2.x << ' ' << model2.y;
//
//    std::cout << model.x;
    E obj(2);
    D &ref = obj;
    ref.printsth();
}