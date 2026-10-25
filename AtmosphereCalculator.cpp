#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;

struct AtmChar {
    double h;
    double H;
    double T;
    double P;
    double p;
    double g;
    double a;
};

AtmChar CalcAtm(double hgeom) {
    AtmChar result;

    result.h = hgeom;

    const double r = 6356767;
    const double g0 = 9.80665;
    const double R = 287.053;
    const double x = 1.4;

    double H0 = 0, T0 = 0, P0 = 0, b = 0, Pl = 0;

    // Геопотенциальная высота
    result.H = r * hgeom / (r + hgeom);

    // Ускорение свободного падения
    result.g = g0 * pow(r / (r + hgeom), 2);

    // Выбор атмосферного слоя

    if (result.H < 0) {
        H0 = -2000;
        T0 = 301.15;
        b = -0.0065;
        P0 = 127774;
    }

    if (result.H >= 0 && result.H < 11000) {
        H0 = 0;
        T0 = 288.15;
        b = -0.0065;
        P0 = 101325;
    }

    if (result.H >= 11000 && result.H < 20000) {
        H0 = 11000;
        T0 = 216.65;
        b = 0;
        P0 = 22632.0;
    }

    if (result.H >= 20000 && result.H < 32000) {
        H0 = 20000;
        T0 = 216.65;
        b = 0.0010;
        P0 = 5474.87;
    }

    if (result.H >= 32000 && result.H < 47000) {
        H0 = 32000;
        T0 = 228.65;
        b = 0.0028;
        P0 = 868.014;
    }

    if (result.H >= 47000 && result.H < 51000) {
        H0 = 47000;
        T0 = 270.65;
        b = 0;
        P0 = 110.906;
    }

    if (result.H >= 51000 && result.H < 71000) {
        H0 = 51000;
        T0 = 270.65;
        b = -0.0028;
        P0 = 66.9384;
    }

    if (result.H >= 71000 && result.H < 85000) {
        H0 = 71000;
        T0 = 214.65;
        b = -0.002;
        P0 = 3.956390;
    }

    if (result.H >= 85000 && result.H < 94000) {
        H0 = 85000;
        T0 = 186.65;
        b = 0;
        P0 = 0.363410;
    }

    // Температура
    result.T = T0 + b * (result.H - H0);

    // Давление
    if (b != 0) {
        Pl = log10(P0) - g0 / (b * R) * log10(result.T / T0);
        result.P = pow(10, Pl);
    }
    else {
        Pl = log10(P0) -
            (0.434294 / (R * result.T)) *
            g0 * (result.H - H0);

        result.P = pow(10, Pl);
    }

    // Плотность воздуха
    result.p = result.P / (R * result.T);

    // Скорость звука
    result.a = sqrt(x * R * result.T);

    return result;
}

string ToScienceA(double value, int precision = 5) {
    if (value == 0.0) return "0.00000";

    stringstream ss;

    int exp = 0;
    double man = value;

    if (fabs(value) >= 1.0) {
        while (fabs(man) >= 10.0) {
            man /= 10.0;
            exp++;
        }
    }
    else {
        while (fabs(man) < 1.0 && man != 0.0) {
            man *= 10.0;
            exp--;
        }
    }

    ss << fixed << setprecision(precision) << man;

    if (exp != 0) {
        ss << " x 10^" << exp;
    }

    return ss.str();
}

string ToScienceB(double value) {
    if (value == 0.0) return "0.000";

    stringstream ss;

    int exp = 0;
    double man = value;

    if (fabs(value) >= 1000.0) {
        while (fabs(man) >= 1000.0) {
            man /= 10.0;
            exp++;
        }
    }
    else if (fabs(value) < 1.0) {
        while (fabs(man) < 100.0 && man != 0.0) {

            if (fabs(man) < 1.0) {
                man *= 1000.0;
                exp -= 3;
            }
            else if (fabs(man) < 10.0) {
                man *= 100.0;
                exp -= 2;
            }
            else if (fabs(man) < 100.0) {
                man *= 10.0;
                exp -= 1;
            }
        }
    }

    ss << fixed << setprecision(3) << man;

    if (exp != 0) {
        ss << " x 10^" << exp;
    }

    return ss.str();
}

void PrintChar(const AtmChar& atm) {
    cout << fixed << setprecision(0)
        << "H = " << atm.H << " (м)" << endl;

    cout << "T = "
        << ToScienceB(atm.T)
        << " (К)" << endl;

    cout << "P = "
        << ToScienceA(atm.P)
        << " (Па)" << endl;

    cout << "p = "
        << ToScienceA(atm.p)
        << " (кг/м^3)" << endl;

    cout << "g = "
        << ToScienceA(atm.g)
        << " (м/с^2)" << endl;

    cout << "a = "
        << ToScienceB(atm.a)
        << " (м/c)" << endl;
}

//  h (м) - Геометрическая высота
//  H (м) - Геопотенциальная высота
//  T (К) - Кинетическая температура
//  P (Па) - Давление воздуха
//  p (кг/м^3) - Плотность воздуха
//  a (м/c) - Скорость звука
//  g (м/с^2) - Ускорение свободного падения

int main() {
    setlocale(LC_ALL, "Rus");

    double h;

    cout << "Input altitude: ";
    cin >> h;

    while (h < -2000 || h > 94000) {
        cout << "Altitude should be [-2000;94000] !" << endl;

        cout << "Input altitude: ";
        cin >> h;
    }

    AtmChar atm = CalcAtm(h);

    cout << "Counted characteristics for "
        << h << " altitude" << endl;

    PrintChar(atm);

    return 0;
}
