#include <cassert>

class Base{
public:
    struct BaseStruct {
        int a;
    };
private:
    const BaseStruct bs_;
public:
    Base(const BaseStruct &s) : bs_(s) {}

    BaseStruct getStruct(){
        return bs_;
    }
};

int main() {
    Base::BaseStruct x = {1};
    Base b(x);
    Base::BaseStruct y = b.getStruct();

    assert(y.a == 1);
}