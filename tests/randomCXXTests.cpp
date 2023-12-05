//#include <iostream>
////#include <omp.h>
//
//using namespace std;
//
//class NrComplex {
//private:
//    double parteReala;
//    double parteImaginara;
//public:
//    NrComplex() {
//        this -> parteReala = 0;
//        this -> parteImaginara = 0;
//    }
//
//    NrComplex(double r, double i) {
//        this -> parteReala = r;
//        this -> parteImaginara = i;
//    }
//
//    void afisare() {
//        cout << this -> parteReala << " + " << this -> parteImaginara << 'i' << '\n';
//    }
//
////    NrComplex operator+(const NrComplex& al_doilea_termen_al_adunarii) {
////        NrComplex rez;
////        rez.parteReala = this -> parteReala + al_doilea_termen_al_adunarii.parteReala;
////        rez.parteImaginara = this -> parteImaginara + al_doilea_termen_al_adunarii.parteImaginara;
////
////        return rez;
////    }
////    void operator+(const NrComplex& al_doilea_termen_al_adunarii){
////        cout << "STO" << '\n';
////    }
//
//    bool operator+(const NrComplex& al_doilea_termen) {
//        cout << "Operator gherla";
//        return false;
//    }
//
//    NrComplex operator-(const NrComplex& al_doilea_termen_al_scaderii){
//        NrComplex rez;
//        rez.parteReala= this -> parteReala - al_doilea_termen_al_scaderii.parteReala;
//        rez.parteImaginara= this -> parteImaginara - al_doilea_termen_al_scaderii.parteImaginara;
//         return rez;
//    }
//
//    bool operator==(const NrComplex& al_doilea_termen_al_comparatiei) {
//        bool rez;
//        int ok=0;
//        if (this->parteReala == al_doilea_termen_al_comparatiei.parteReala && this->parteImaginara== al_doilea_termen_al_comparatiei.parteImaginara)
//            ok = 1;
//
//        if (ok == 0) {
//            return false;
//        } else {
//            return true;
//        }
//    }
//
//    NrComplex operator+(double numar_real) {
//        NrComplex rez;
//        rez.parteReala = this -> parteReala + numar_real;
//        rez.parteImaginara = this -> parteImaginara;
//
//        return rez;
//    }
//
//    double getParteReala() {
//        return this -> parteReala;
//    }
//
//    double getParteImaginara() {
//        return this -> parteImaginara;
//    }
//
//    void setParteReala(double p) {
//        this -> parteReala = p;
//    }
//
//    void setParteImaginara(double i) {
//        this -> parteImaginara = i;
//    }
//};
//
//// OPERATOR+ scris in afara clasei
//NrComplex operator+(double numar_real, NrComplex& numar) {
//    NrComplex rez(numar_real + numar.getParteReala(), numar.getParteImaginara());
//    return rez;
//
//    // SAU
////     NrComplex rez;
////     rez.setParteReala(numar_real + numar.getParteReala());
////     rez.setParteImaginara(numar.getParteImaginara());
////     return rez;
//}
//
//int main() {
//
//    NrComplex z1(3, 4), z2(-1, 1);
////    NrComplex suma = z1 + z2; //z1.adunare(z2)
////    suma.afisare();
////    NrComplex scadere= z1-z2;
////    scadere.afisare();
////    NrComplex z3(3, 4);
////    if (z1 == z3) {
////        cout << "Egale";
////    }
////    else {
////        cout << "Nu sunt egale";
////    }
////    int x = 4;
////    NrComplex z5 = z1;
////    z5.afisare();
////    z1 + z2; // z1.adunare(z2)
//
//    // z[1] partea reala
//    // z[2] partea imaginara
////    if (z1 + z2) {
////        ;
////    }
////    else {
////        cout << '\n' << "Am intrat pe else";
////    }
//    NrComplex z6 = z1 + 4;
//    z6.afisare();
//
//}
//
//
//
//#include <iostream>
//#include<string>
//#include <cmath>
//
//using namespace std;
//
//class Vector {
//private:
//    const int id;
//    int n;
//    float* v;
//    static int nr;
//    float modul;
//
//public:
//
//    Vector() : id(++nr) {
//        this-> n = 0;
//        this->modul = 0.0f;
//        this->v = nullptr;
//
//
//    }
//
//    // constructor cu un singur parametru
//    Vector(int n) : id(++nr) {
//        this -> n = n;
//        this -> v = new float[this -> n];
//
//        for (int i = 0; i < this -> n; i++) {
//            this -> v[i] = 0;
//        }
//        this -> modul = 0;
//    }
//
//    // OBS 1
//    // nu ai facut nimic cu parametrul modul
//    // putea fi sters
//    Vector(int n, float modul, float* v2) : id(++nr) {
//        this->n = n;
//        int suma = 0;
//
//        this->v = new float[this->n];
//        for (int i = 0; i < this->n; i++) {
//            this->v[i] = v2[i];
//        }
//        for (int i = 0;i < this->n; i++)
//            suma = suma + v[i] * v[i];
//        this->modul = sqrt(suma);
//        delete[] v2;
//    }
//
//    Vector(const Vector& copie) : id(copie.id) {
//        this->n = copie.n;
//        this->modul = copie.modul;
////        if (copie.v != nullptr) {
////            this->v = new float[this->n];
////            memcpy(this->v, copie.v, this->n);
////        }
////        else {
////            this->v = nullptr;
////        }
//
//        // Nu face cu memcpy
//        // fa ca la seminar
//        if (copie.v != nullptr) {
//            this -> v = new float[this -> n];
//            for (int i = 0; i < this -> n; i++) {
//                this -> v[i] = copie.v[i];
//            }
//        }
//        else {
//            this -> v = nullptr;
//        }
//    }
//
//    // OPERATOR=
//    // E important
//    Vector& operator=(const Vector& copie) {
//        if (this == &copie) {
//            return *this;
//        }
//
//        Vector rez;
//        this->n = copie.n;
//        this->modul = copie.modul;
//
//        if (copie.v != nullptr) {
//            this -> v = new float[this -> n];
//            for (int i = 0; i < this -> n; i++) {
//                this -> v[i] = copie.v[i];
//            }
//        }
//        else {
//            this -> v = nullptr;
//        }
//
//        return *this;
//    }
//
//    void afisare() {
//        cout << endl << "vectorul " << "de " << this->n << " elemente " << " cu urmatoarele valori: " << endl;
//        for (int i = 0; i < this->n; i++) {
//            cout << this->v[i] << " ";
//        }
//        cout << "si modulul egal cu " << this->modul;
//    }
//
//    ~Vector() {
//        delete this->v;
//    }
//
//    // OBS 5
//    // Adaugam si modulul lui rez
//    Vector operator+(const Vector& v3) {
//        Vector rez;
//        rez.n = n;
//        rez.v = new float[rez.n];
//        for (int i = 0;i < this->n; i++) {
//            rez.v[i] = this->v[i] + v3.v[i];
//        }
//
//        // sqrt e din #include <cmath>
//        rez.modul = sqrt((this -> modul) * (this -> modul) + v3.modul * v3.modul);
//
//        return rez;
//    }
//
//
//    Vector operator-(const Vector & v3) {
//
//        Vector rez;
//        rez.n = n;
//        rez.v = new float[rez.n];
//        for (int i = 0;i < this->n; i++) {
//            rez.v[i] = this->v[i] - v3.v[i];
//        }
//        return rez;
//
//    }
//
//    // OBS 3
//    // ar trebui float multiplicator ca inmultesti cu floaturi
//
//    // IMPORTANT
//    // Desi e corect, aici ar fin mai bine sa returnez un vector
//    // Altfel operatorul ar trebui apelat la modul obj * 7;
//    // As vrea obj = obj * 7, care nu merge daca returnez void
////    void operator*=(int multiplicator) {
////        for (int i = 0; i < this->n; i++)
////            this->v[i] *= multiplicator;
////
////    }
//    Vector operator*=(float multiplicator) {
//        Vector rez(this -> n);
//        for (int i = 0; i < this -> n; i++){
//            rez.v[i] = this -> v[i] * multiplicator;
//        }
//
//        rez.modul *= multiplicator;
//        return rez;
//    }
//
//    // Adaug si operator * folosind *=
//    Vector operator*(float multiplicator) {
//        return *this *= multiplicator;
//    }
//
//    // OBS 6
//    // Aici nu am de ce sa returnez adresa
//    //float& operator[](int index) {
//    float operator[](int index) {
//        if (index >= 0 && index <this->n) {
//            return this->v[index];
//        }
//        else {
//            throw "Indexul nu se afla in interval";
//        }
//    }
//
//
//    bool operator>(const Vector& v2)
//    {
//        return this->modul > v2.modul;
//    }
//    bool operator<(const Vector& v2) {
//
//        return this->modul < v2.modul;
//    }
//    bool operator==(const Vector& v2) {
//        return this->modul == v2.modul;
//    }
//    static int getnr() {
//        return nr;
//    }
//    float getModul() {
//        return this->modul;
//    }
//
//
//};
//
//int Vector:: nr = 0;
//
//// OBS 4
//// Fa int main() ca pe unele compilatoare nu merge cu void
//int main() {
//    Vector obiect1(3, 50, new float[3] {1, 2, 3});
//    Vector obiect3(3, 50, new float[3] {-4, -5, -6});
//
//    Vector rezultat = obiect1 + obiect3;
//    rezultat.afisare();
//    Vector* obiect2;
//    Vector rezultat2 = obiect1 - obiect3;
//    rezultat2.afisare();
//    obiect1 *= 2;
//    obiect1.afisare();
//    cout << endl;
//    if (obiect1 > obiect3) {
//        cout << 1;
//    }
//    else cout << 0;
//
//
//    // SCRIS DE MINE
//
//    // test copy constructor
//    Vector obj3(obiect1);
//    obj3.afisare();
//
//    obj3 = obiect1 * 3;
//    obj3.afisare();
//
//    // problema pe care o ziceam eu
//    int dim = 3;
//    Vector* vectoriDeBaza = new Vector[dim];
//    for (int i = 0; i < dim; i++) {
//        float* elemente = new float[dim];
//        for (int j = 0; j < dim; j++) {
//            if (i != j) {
//                elemente[j] = 0;
//            }
//            else if (i == j) {
//                elemente[j] = 1;
//            }
//        }
//
//        Vector e_i(dim, -1, elemente);
//        vectoriDeBaza[i] = e_i;
//    }
//
//    Vector suma(dim);
//    for (int i = 0; i < dim; i++) {
//        vectoriDeBaza[i].afisare();
//        suma = suma + vectoriDeBaza[i] * obiect1[i];
//    }
//
//    obiect1.afisare();
//    suma.afisare();
//
//
//    // OBS 7
//    // Sa nu uiti sa dezaloci memoria la final de main
//    return 0;
//}

#include <iostream>
#include <omp.h>

int main() {
    int tid, nthreads;

#pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        std::cout << "Thread " << tid << " of " << nthreads << '\n';
    }

    return 0;
}