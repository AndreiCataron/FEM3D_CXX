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

//#include <iostream>
//#include <omp.h>
//
//int main() {
//    int tid, nthreads;
//
//#pragma omp parallel private(tid)
//    {
//        tid = omp_get_thread_num();
//        nthreads = omp_get_num_threads();
//
//        std::cout << "Thread " << tid << " of " << nthreads << '\n';
//    }
//
//    return 0;
//}

//#include<iostream>
//#include<string>
//
//using namespace std;
//
//class festivale {
//private:
//    const int id_festival;
//
//public:
//    string numeFestival;
//    string main_artist;
//    float *pret;
//
//public:
//
//    festivale() : id_festival(0), numeFestival("Electric Castle"), main_artist("Necunoscut"), pret(nullptr) {}
//
//
//    festivale(int id, string numeFestival, string main_artist, float *pret) : id_festival(id),
//                                                                              numeFestival(numeFestival),
//                                                                              main_artist(main_artist), pret(pret) {}
//
//    festivale(const festivale& altFest) : id_festival(altFest.id_festival) {
//        numeFestival = altFest.numeFestival;
//        main_artist = altFest.main_artist;
//
//        pret = nullptr;
//    }
//
//    ~festivale() {
//        cout << endl << "A fost apelat destructorul clasei de baza." << endl;
//        delete pret;
//    }
//
//
//    float* getPret() {
//        return pret;
//    }
//
//};
//
//void afisareCampuri(festivale *obiect) {
//    cout << endl << "Festivalul " << obiect->numeFestival << " cu main artistul " << obiect->main_artist;
//}
//
//enum TipMuzica { Simfonica=1, Camerata=2, Corala=3 };
//
//class FestivalMuzicaClasica : public festivale {
//private:
//    int* locuriPrimulRand;
//    char* numeOrchestra;
//    const bool esteInstrumental;
//    TipMuzica tipMuzica;
//    int durataConcert;
//    int nrInstrumente;
//    int nrBileteVandute;
//    int nrLocuriPrimulRand;
//    float pretBilet;
//
//public:
//    FestivalMuzicaClasica() : festivale(), locuriPrimulRand(nullptr), esteInstrumental(true), tipMuzica(Simfonica), durataConcert(120), nrInstrumente(50), nrBileteVandute(0), nrLocuriPrimulRand(10), pretBilet(50.0f) {
//        this -> numeOrchestra = new char[100];
//        strcpy(this -> numeOrchestra, "Orchestra Necunoscuta");
//        this -> locuriPrimulRand = new int[nrLocuriPrimulRand];
//        for (int i = 0; i < nrLocuriPrimulRand; i++) {
//            this -> locuriPrimulRand[i] = i + 1;
//        }
//    }
//
//    FestivalMuzicaClasica(int id, string numeFestival, string main_artist, float *pret, int* locuri, char* nume, bool instrumental, TipMuzica tip, int durata, int nrInst, int nrBilete, int nrScaune, float pret2) :
//                          festivale(id, numeFestival, main_artist, pret), esteInstrumental(instrumental), tipMuzica(tip), durataConcert(durata), nrInstrumente(nrInst), nrBileteVandute(nrBilete), nrLocuriPrimulRand(nrScaune), pretBilet(pret2) {
//        this -> locuriPrimulRand = new int[this -> nrLocuriPrimulRand];
//        for (int i = 0; i < this -> nrLocuriPrimulRand; i++) {
//            this -> locuriPrimulRand[i] = locuri[i];
//        }
//
//        this -> numeOrchestra = new char[strlen(nume) + 1];
//        memcpy(this -> numeOrchestra, nume, strlen(nume) + 1);
//
//        delete[] locuri;
//    }
//
//    FestivalMuzicaClasica(TipMuzica tip, int durata, int nrInst, float pret2) :
//            festivale(), esteInstrumental(false), tipMuzica(tip), durataConcert(durata), nrInstrumente(nrInst), nrBileteVandute(0), nrLocuriPrimulRand(0), pretBilet(pret2) {
//        this -> locuriPrimulRand = nullptr;
//
//        this -> numeOrchestra = nullptr;
//    }
//
//    FestivalMuzicaClasica(const FestivalMuzicaClasica& altFestival) : festivale(altFestival), esteInstrumental(altFestival.esteInstrumental), tipMuzica(altFestival.tipMuzica), durataConcert(altFestival.durataConcert), nrInstrumente(altFestival.nrInstrumente), nrBileteVandute(altFestival.nrBileteVandute), nrLocuriPrimulRand(altFestival.nrLocuriPrimulRand), pretBilet(altFestival.pretBilet) {
//        this -> locuriPrimulRand = new int[altFestival.nrLocuriPrimulRand];
//        for (int i = 0; i < this -> nrLocuriPrimulRand; i++) {
//            this -> locuriPrimulRand[i] = altFestival.locuriPrimulRand[i];
//        }
//
//        this -> numeOrchestra = new char[strlen(altFestival.numeOrchestra) + 1];
//        strcpy(this -> numeOrchestra, altFestival.numeOrchestra);
//    }
//
//    ~FestivalMuzicaClasica() {
//        cout << "A fost apelat destructorul derivatei." << endl;
//        delete[] this -> locuriPrimulRand;
//        delete[] this -> numeOrchestra;
//    }
//
//    FestivalMuzicaClasica& operator=(const FestivalMuzicaClasica& altFestival) {
//        if (this != &altFestival) {
//            delete[] this -> locuriPrimulRand;
//            delete[] this -> numeOrchestra;
//
//            this -> locuriPrimulRand = new int[altFestival.nrLocuriPrimulRand];
//            for (int i = 0; i < altFestival.nrLocuriPrimulRand; i++) {
//                this -> locuriPrimulRand[i] = altFestival.locuriPrimulRand[i];
//            }
//
//            numeOrchestra = new char[strlen(altFestival.numeOrchestra) + 1];
//            strcpy(numeOrchestra, altFestival.numeOrchestra);
//
//            this -> tipMuzica = altFestival.tipMuzica;
//            this -> durataConcert = altFestival.durataConcert;
//            this -> nrInstrumente = altFestival.nrInstrumente;
//            this -> nrBileteVandute = altFestival.nrBileteVandute;
//            this -> nrLocuriPrimulRand = altFestival.nrLocuriPrimulRand;
//            this -> pretBilet = altFestival.pretBilet;
//
//        }
//        return *this;
//    }
//
//
//
//    char* getNumeOrchestra() {
//        return this -> numeOrchestra;
//    }
//
//    int getPret() {
//        return this -> pretBilet;
//    }
//
//    void setNumeOrchestra(const char* nume) {
//        delete[] this -> numeOrchestra;
//        this -> numeOrchestra = new char[strlen(nume) + 1];
//        strcpy(this -> numeOrchestra, nume);
//    }
//
//    int operator() () {
//        return this -> nrBileteVandute;
//    }
//
//    void operator+=(float scumpire) {
//        this -> pretBilet += scumpire;
//    }
//
//    void afisare() {
//        afisareCampuri(this);
//
//        cout << ". Acesta este un festival de muzica de muzica clasica de tip " << this -> tipMuzica << ". Concertul va dura " << this -> durataConcert << " minute, va avea " <<
//        this -> nrInstrumente << " tipuri de instrumente. In sala sunt disponibile " << this -> nrLocuriPrimulRand << " locuri pe primul rand. Va canta orchestra ";
//        if(numeOrchestra != nullptr){
//            cout << numeOrchestra << ". \n";
//        }
//        else {
//            cout << " fara nume." << endl;
//        }
//    }
//};
//
//int main() {
//
//    float* pret = new float[100.0f];
//    int* locuri = new int[5] {1, 2, 3, 4, 5};
//    char* nume1 = new char[strlen("George Enescu")];
//    strcpy(nume1, "George Enescu");
//
//    FestivalMuzicaClasica festival_default();
//
//    FestivalMuzicaClasica george_enescu(0, nume1, "Orchestra Simfonica Bucuersti", pret, locuri, "OSB", true, Simfonica, 100, 10, 200, 5, 100);
//
//    FestivalMuzicaClasica* copie = new FestivalMuzicaClasica(george_enescu);
//
//    copie -> afisare();
//
//    cout << "S-au vandut " << (*copie)() << " bilete. Canta orchestra " << copie -> getNumeOrchestra() << endl;
//
//    george_enescu.setNumeOrchestra("Filarmonica Constanta");
//    cout << copie -> getNumeOrchestra() << " nu mai poate ajunge la festival. Va fi inlocuita de " << george_enescu.getNumeOrchestra() << endl;
//
//    george_enescu += 50;
//    cout << "S-a scumpit biletul! Vechiul pret era " << copie -> getPret() << ". Noul pret este " << george_enescu.getPret();
//
//    FestivalMuzicaClasica patru_parametrii(Corala, 100, 3, 200);
//    patru_parametrii.afisare();
//
//    return 0;
//}

//#include<iostream>
//#include <cfloat>
//using namespace std;
//
//// Function for Machine Epsilon with an
//// initial value provided as EPS.
//void machineEpsilon(float EPS)
//{
//    // taking a floating type variable
//    float prev_epsilon;
//
//    // run until condition satisfy
//    while ((1+EPS) != 1)
//    {
//        // copying value of epsilon into previous epsilon
//        prev_epsilon = EPS;
//
//        // dividing epsilon by 2
//        EPS /=2;
//    }
//
//    // print output of the program
//    cout << "Machine Epsilon is : " << prev_epsilon << endl;
//}
//#include <iostream>
//
//constexpr int factorial_cxx14(int n)
//{
//    int res = 1;
//    while (n > 1)
//        res *= n--;
//    return res;
//}
//
//int main()
//{
//    int n;
//    std::cin >> n;
//    std::cout << factorial_cxx14(n);
//    return 0;
//}

#include <iostream>

using namespace std;

class Poligon{
private:
    int nr_laturi;
protected:
    double arie;
public:
    Poligon() {
        cout << "Constructor Poligon" << endl;
    }

    virtual void calculeaza_arie() = 0;
    double get_arie() {
        return this -> arie;
    };

    ~Poligon() {
        cout << "Destructor Poligon" << endl;
    }
};

class Triunghi : public Poligon {
private:
    double a, b, c;
public:
    Triunghi(double a, double b, double c) {
        this -> a = a;
        this -> b = b;
        this -> c = c;
        cout << "Constructor Triunghi" << endl;
    }

    void calculeaza_arie() override{
        cout << 2;
    }
};

class Patrulater : public Poligon{
public:
    Patrulater() {
        cout << "Connstructor Patrulater" << endl;
    }

    ~Patrulater() {
        cout << "Destructor Patrulater" << endl;
    }
};

class Patrat : public Patrulater{
private:
    double l;
public:
    Patrat(double latura) : Patrulater(){
        cout << "Constructor Patrat" << endl;
        this -> l = latura;
    }

    void calculeaza_arie() override {
        this -> arie =  l * l;
    }

    ~Patrat() {
        cout << "Destructor Patrat" << endl;
    }
};

//class Desen{
//private:
//    int nr_poligoane;
//    Poligon* poligoane;
//public:
//    Desen(int nrpoligoane, Poligon* vector_poligoane){
//        this -> nr_poligoane = nrpoligoane;
//        this -> poligoane = new Poligon[this -> nr_poligoane];
//        for (int i = 0; i < this -> nr_poligoane; i++) {
//            this -> poligoane[i] = vector_poligoane[i];
//        }
//        delete[] vector_poligoane;
//    }
//};

int main() {
    Patrat p(2);
    Triunghi t(1, 2, 3);
    Poligon* poli1 = &p;
    Poligon* poli2 = &t;

    Poligon* v[3];

    vector<Poligon*> vec;
    vec.push_back(poli1);
    vec.push_back(poli2);

   //vector<int> v = {1, 2, 3};
    cout << v[2] << endl;

}