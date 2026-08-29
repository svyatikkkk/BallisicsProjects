#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

const double PI = 3.141592653589793;
const double MU = 398600.4418; 
const double Rad2Deg = 180.0 / PI;
const double Err = 1e-3; 

static double NormaliseAngles(double angle) {
    while (angle < 0) angle = angle + 360;
    while (angle >= 360) angle = angle - 360;
    return angle;
}

struct KeplerElements {
    double SMA; 
    double ECC; 
    double INC; 
    double RAAN; 
    double AOP; 
    double TA; 
    double EA; 
    double MA; 
    double u; 
    string OrbitType;
    string OrbitShape;
};

struct Vector {
    double x; 
    double y; 
    double z; 
};

int main() {
	setlocale(LC_ALL, "Russian");
	Vector V{};    
    Vector r{};
    Vector c{};
    KeplerElements result{};

    double r_mag, p, a; 
    double V_mag, V_r; 
    double c_mag; 
    double Energy; 
    double e;
    double i, OMG, omg, u, theta, E, M; 
    double sin_i, cos_i, sin_OMG, cos_OMG, sin_u, cos_u, sin_TA = 0, cos_TA = 0;
    double sin_E, cos_E;
    string temp;

    cout << "Программа по определению элементов Кеплера по радиус-вектору и вектору скорости и наоборот" << endl;
    cout << "1 - Определение элементов Кеплера по радиус-вектору и вектору скорости" << endl;
    cout << "2 - Определение радиус-вектора и вектора скорости по элементам Кеплера" << endl;
    cout << endl;
    cout << "Выберите режим работы программы (1 или 2): ";

    int choice;
    cin >> choice;

    switch (choice) {
    case 1: {
        cout << "Вы выбрали режим определения элементов Кеплера по радиус-вектору и вектору скорости" << endl;
        cout << endl;

        ifstream input_file("input.txt");
        input_file >> temp >> V.x;
        input_file >> temp >> V.y;
        input_file >> temp >> V.z;
        input_file >> temp >> r.x;
        input_file >> temp >> r.y;
        input_file >> temp >> r.z;
        input_file.close();

        r_mag = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);

        V_mag = sqrt(V.x * V.x + V.y * V.y + V.z * V.z);

        c.x = r.y * V.z - r.z * V.y;
        c.y = r.z * V.x - r.x * V.z;
        c.z = r.x * V.y - r.y * V.x;
        c_mag = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);

        p = (c_mag * c_mag) / MU;

        Energy = (V_mag * V_mag) / 2.0 - MU / r_mag;
        a = -MU / (2.0 * Energy);
        result.SMA = a;

        e = sqrt(1.0 - p / a);
        if (fabs(e) < Err) {
            result.OrbitShape = "Circular";
        }
        else if (fabs(e - 1.0) < Err) {
            result.OrbitShape = "Parabolical";
        }
        else if ((e > Err) && (e < 1.0)) {
            result.OrbitShape = "Elliptical";
        }
        else result.OrbitShape = "Hyperbolical";

        result.ECC = e;

        cos_i = c.z / c_mag;
        if (cos_i < -1) cos_i = -1;
        if (cos_i > 1) cos_i = 1;
        sin_i = sqrt(1 - cos_i * cos_i);
        i = acos(cos_i);
        result.INC = Rad2Deg * i;

        if (fabs(cos_i - 1) < Err) result.OrbitType = "Equatorial";
        else if (fabs(cos_i + 1) < Err) result.OrbitType = "Equatorial_reverse";
        else if (fabs(cos_i) < Err) result.OrbitType = "Polar";
        else result.OrbitType = "Inclined";

        if (sin_i > Err) {
            sin_OMG = c.x / (c_mag * sin_i);
            cos_OMG = -c.y / (c_mag * sin_i);
            OMG = atan2(sin_OMG, cos_OMG);
            result.RAAN = OMG * Rad2Deg;
        }
        else {
            OMG = 0.0;
            result.RAAN = OMG * Rad2Deg;
        }

        if (sin_i > Err) {
            sin_u = r.z / (r_mag * sin_i);
            cos_u = (r.x * cos_OMG + r.y * sin_OMG) / r_mag;
            u = atan2(sin_u, cos_u);
            result.u = u * Rad2Deg;
        }
        else {
            u = atan2(r.y, r.x);
            result.u = u * Rad2Deg;
        }

        if (fabs(e) > Err) {
            V_r = (r.x * V.x + r.y * V.y + r.z * V.z) / r_mag;
            sin_TA = r_mag * V_r / (e * sqrt(MU * p));
            cos_TA = (p - r_mag) / (e * r_mag);
            theta = atan2(sin_TA, cos_TA);
            result.TA = theta * Rad2Deg;
        }
        else {
            theta = 0.0;
            result.TA = theta * Rad2Deg;
        }

        if (fabs(e) > Err) {
            sin_E = (sqrt(1 - e * e) * sin_TA) / (1 + e * cos_TA);
            cos_E = (e + cos_TA) / (1 + e * cos_TA);
            E = atan2(sin_E, cos_E);
            result.EA = E * Rad2Deg;
        }
        else {
            E = theta;
            result.EA = E * Rad2Deg;
        }

        if (fabs(e) > Err) {
            M = E - e * sin(E);
            result.MA = M * Rad2Deg;
        }
        else {
            M = 0;
            result.MA = M * Rad2Deg;
        }

        if (fabs(e) > Err) {
            omg = u - theta;
            result.AOP = omg * Rad2Deg;
        }
        else {
            omg = u;
            result.AOP = omg * Rad2Deg;
        }

        cout << "Вывести результаты в консоль или в файл output.txt? (1 - консоль, 2 - .txt): ";
        int choice2;
        cin >> choice2;

        switch (choice2) {
        case 1:
            cout << "KeplerElements" << endl;
            cout << "SMA        = " << result.SMA << " km" << endl;
            cout << "ECC        = " << result.ECC << endl;
            cout << "INC        = " << NormaliseAngles(result.INC) << " deg" << endl;
            cout << "RAAN       = " << NormaliseAngles(result.RAAN) << " deg" << endl;
            cout << "AOP        = " << NormaliseAngles(result.AOP) << " deg" << endl;
            cout << "TA         = " << NormaliseAngles(result.TA) << " deg" << endl;
            cout << "EA         = " << NormaliseAngles(result.EA) << " deg" << endl;
            cout << "MA         = " << NormaliseAngles(result.MA) << " deg" << endl;
            cout << "u          = " << NormaliseAngles(result.u) << " deg" << endl;
            cout << "Orbit type = " << result.OrbitType << endl;
            cout << "Orbit shape = " << result.OrbitShape << endl;
            break;

        case 2: {
            ofstream output_file("output.txt");
            output_file << "SMA " << result.SMA << endl;
            output_file << "ECC " << result.ECC << endl;
            output_file << "INC " << NormaliseAngles(result.INC) << endl;
            output_file << "RAAN " << NormaliseAngles(result.RAAN) << endl;
            output_file << "AOP " << NormaliseAngles(result.AOP) << endl;
            output_file << "TA " << NormaliseAngles(result.TA) << endl;
            output_file << "EA " << NormaliseAngles(result.EA) << endl;
            output_file << "MA " << NormaliseAngles(result.MA) << endl;
            output_file << "u " << NormaliseAngles(result.u) << endl;
            output_file << "OrbitType " << result.OrbitType << endl;
            output_file << "OrbitShape " << result.OrbitShape << endl;
            output_file.close();
            cout << "Результаты записаны в файл output.txt" << endl;
            break;
        }

        default:
            cout << "Неверный выбор! Введите 1 или 2." << endl;
            break;
        }

        break;
    }

    case 2:
        cout << "Вы выбрали режим определения радиус-вектора и вектора скорости по элементам Кеплера" << endl;
        cout << endl;
        cout << "Режим 2 пока не реализован." << endl;
        break;

    default:
        cout << "Неверный выбор! Введите 1 или 2." << endl;
        break;
    }

    return 0;
}