#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

double V, psi, x, z, t;

double a_atm, ro_atm, p_atm, g_atm, T_atm;

int nza = 0;
double h = 17000;
double v0 = 440;
double t0 = 0;
double psi0 = 0;
double x0 = 0;
double z0 = 0;
double P0 = 25543;
double Sm = 0.388;
double Sa = 0.13;
double poN = 101325;
double t_ayt = 3.59;
double PI = 3.14159265358979323846;
double k = 2.24;

void CalcAtm(double h) {
    a_atm = 295.070;
    ro_atm = 0.142301;
    p_atm = 8849.72;
    g_atm = 9.7545;
    T_atm = 216.650;
}

double Interpolate_cR(double a) {
    double MachVals[11] = { 1.4, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5 };
    double cR[11] = { 0.0391, 0.0345, 0.0306, 0.0276, 0.0254, 0.0236, 0.0224, 0.0216, 0.0214, 0.0210 , 0.0207 };
    double cRa = 0.01;
    double M = a;
    for (int i = 0; i < 10; i++) {
        if (MachVals[i] <= M && M < MachVals[i + 1]) {
            cRa = cR[i] + ((M - MachVals[i]) * (cR[i + 1] - cR[i]) / (MachVals[i + 1] - MachVals[i]));
        }
    }
    return cRa;
}

double Interpolate_Ca(double alf, double M) {
    double AlphaVals[19] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20 };
    double ca_1[19] = { 0.0462, 0.0478, 0.0500, 0.0525, 0.0573, 0.0620, 0.0697, 0.0770, 0.0867, 0.0956, 0.1069, 0.1183, 0.1268, 0.1397, 0.1533, 0.1678, 0.1830, 0.2145, 0.2367 };
    double ca_5[19] = { 0.0310, 0.0320, 0.0330, 0.0350, 0.0380, 0.0410, 0.0440, 0.0490, 0.0540, 0.0600, 0.0660, 0.0745, 0.0819, 0.0913, 0.1018, 0.1131, 0.1251, 0.1540, 0.1835 };

    double Cxa = 0;
    double abs_alf = abs(alf);

    if (abs_alf <= AlphaVals[0]) {
        double c_14 = ca_1[0];
        double c_5 = ca_5[0];
        Cxa = (c_14 + (M - 1.4) * (c_5 - c_14) / (5 - 1.4));
    }
    else if (abs_alf >= AlphaVals[18]) {
        double c_14 = ca_1[18];
        double c_5 = ca_5[18];
        Cxa = (c_14 + (M - 1.4) * (c_5 - c_14) / (5 - 1.4));
    }
    else {
        for (int i = 0; i < 18; i++) {
            if (AlphaVals[i] <= abs_alf && abs_alf < AlphaVals[i + 1]) {
                double c_14 = ca_1[i] + (abs_alf - AlphaVals[i]) * (ca_1[i + 1] - ca_1[i]) / (AlphaVals[i + 1] - AlphaVals[i]);
                double c_5 = ca_5[i] + (abs_alf - AlphaVals[i]) * (ca_5[i + 1] - ca_5[i]) / (AlphaVals[i + 1] - AlphaVals[i]);
                Cxa = (c_14 + (M - 1.4) * (c_5 - c_14) / (5 - 1.4));
                break;
            }
        }
    }

    return Cxa;
}

double Mach(double V) {
    return V / a_atm;
}

double nRalf(double P, double cR, double M, double A, double m) {
    return (P * (PI / 180) + cR * M * M * A) / (m * g_atm);
}

double CalcAlpha(int nza, double nR) {
    if (nza >= 0) {
        return (nza + 2) / nR;
    }
    else {
        return -((-(double)nza) + 2) / nR;
    }
}

double dV_dt(double m, double M, double Cxa, double P, double A) {
    return (1 / m) * (P - Cxa * M * M * A);
}

double dpsi_dt(int nza, double V) {
    return (-g_atm * nza / V);
}

double dx_dt(double V, double psi) {
    return V * cos(psi);
}

double dz_dt(double V, double psi) {
    return (-V * sin(psi));
}

double CalcMass(double t) {
    const double m0 = 129;
    const double dm = 17.7;
    if (t < 0) return m0;
    if (t > t_ayt) return CalcMass(t_ayt);
    return m0 - dm * t;
}

void CalcDerivatives(double P, double A, double m,
    double& dV, double& dpsi,
    double& dx, double& dz) {

    double M_val = Mach(V);
    double cR = Interpolate_cR(M_val);
    double nR = nRalf(P, cR, M_val, A, m);
    double alf = CalcAlpha(nza, nR);
    double Cxa = Interpolate_Ca(alf, M_val);

    dV = dV_dt(m, M_val, Cxa, P, A);
    dpsi = dpsi_dt(nza, V);
    dx = dx_dt(V, psi);
    dz = dz_dt(V, psi);
}

void RungeKutta(double dt, double P, double A) {

    double m_current = CalcMass(t);

    double V0 = V, psi0_local = psi, x0_local = x, z0_local = z, t0_local = t;

    double dV1, dpsi1, dx1, dz1;
    CalcDerivatives(P, A, m_current, dV1, dpsi1, dx1, dz1);

    V = V0 + 0.5 * dt * dV1;
    psi = psi0_local + 0.5 * dt * dpsi1;
    x = x0_local + 0.5 * dt * dx1;
    z = z0_local + 0.5 * dt * dz1;
    t = t0_local + 0.5 * dt;
    double m_mid1 = CalcMass(t);

    double dV2, dpsi2, dx2, dz2;
    CalcDerivatives(P, A, m_mid1, dV2, dpsi2, dx2, dz2);

    V = V0 + 0.5 * dt * dV2;
    psi = psi0_local + 0.5 * dt * dpsi2;
    x = x0_local + 0.5 * dt * dx2;
    z = z0_local + 0.5 * dt * dz2;
    t = t0_local + 0.5 * dt;
    double m_mid2 = CalcMass(t);

    double dV3, dpsi3, dx3, dz3;
    CalcDerivatives(P, A, m_mid2, dV3, dpsi3, dx3, dz3);

    V = V0 + dt * dV3;
    psi = psi0_local + dt * dpsi3;
    x = x0_local + dt * dx3;
    z = z0_local + dt * dz3;
    t = t0_local + dt;
    double m_next = CalcMass(t);

    double dV4, dpsi4, dx4, dz4;
    CalcDerivatives(P, A, m_next, dV4, dpsi4, dx4, dz4);

    V = V0 + (dt / 6.0) * (dV1 + 2 * dV2 + 2 * dV3 + dV4);
    psi = psi0_local + (dt / 6.0) * (dpsi1 + 2 * dpsi2 + 2 * dpsi3 + dpsi4);
    x = x0_local + (dt / 6.0) * (dx1 + 2 * dx2 + 2 * dx3 + dx4);
    z = z0_local + (dt / 6.0) * (dz1 + 2 * dz2 + 2 * dz3 + dz4);
    t = t0_local + dt;
}

int main() {
    setlocale(LC_ALL, "Russian");

    string filename1 = "nza = " + to_string(nza) + "_RK4.txt";
    string filename2 = "nza = " + to_string(nza) + " для отчета_RK4.txt";

    ofstream Predv, Results;
    Predv.open(filename1);
    Results.open(filename2);

    cout << "ИНТЕГРИРОВАНИЕ МЕТОДОМ РУНГЕ-КУТТА 4" << endl;
    cout << "Начальные данные:" << endl;
    cout << "Высота: " << h << " м" << endl;
    cout << "Начальная скорость: " << fixed << setprecision(1) << v0 << " м/с" << endl;
    cout << "Начальный угол: " << fixed << setprecision(1) << psi0 << " рад" << endl;
    cout << "Перегрузка nza: " << nza << endl;
    cout << "Время работы двигателя t_ayt: " << fixed << setprecision(2) << t_ayt << " с" << endl;
    cout << "-------------------------------------------------------------------------------------" << endl;
    cout << setw(7) << "t"
        << setw(8) << "m"
        << setw(10) << "M"
        << setw(11) << "V"
        << setw(13) << "psi"
        << setw(10) << "x"
        << setw(12) << "z"
        << setw(11) << "alf" << endl;
    cout << "-------------------------------------------------------------------------------------" << endl;
    Predv << "t m M alf cR cA V psi x z" << endl;
    Results << "t m M alf cR cA V psi x z" << endl;

    CalcAtm(h);
    double A = (ro_atm * a_atm * a_atm * Sm / 2);
    double P = P0 + Sa * (poN - p_atm);
    double Vmin = k * a_atm;
    double Mmin = Vmin / a_atm;

    V = v0;
    psi = psi0;
    x = x0;
    z = z0;
    t = 0;

    double m_val = CalcMass(t);
    double M_val = Mach(V);
    double cR = Interpolate_cR(M_val);
    double nR = nRalf(P, cR, M_val, A, m_val);
    double alf = CalcAlpha(nza, nR);
    double Cxa = Interpolate_Ca(alf, M_val);

    Predv << t << " "
        << m_val << " "
        << M_val << " "
        << alf << " "
        << cR << " "
        << Cxa << " "
        << V << " "
        << psi * 180 / PI << " "
        << x << " "
        << z << endl;

    Results << t << " "
        << m_val << " "
        << M_val << " "
        << alf << " "
        << cR << " "
        << Cxa << " "
        << V << " "
        << psi * 180 / PI << " "
        << x << " "
        << z << endl;

    cout << fixed << setprecision(1)
        << setw(8) << t
        << fixed << setprecision(4)
        << setw(10) << m_val
        << setw(10) << M_val
        << setw(12) << V
        << setw(10) << psi * 180 / PI
        << setw(12) << x
        << setw(12) << z
        << setw(10) << alf << endl;

    bool t_ayt_reached = false;
    double dt;
    double prev_V, prev_psi, prev_x, prev_z, prev_t;

    for (int i = 0; i < 44; i++) {
        if (i < 35) {
            dt = 0.1;
        }
        else {
            dt = 0.01;
        }

        prev_V = V;
        prev_psi = psi;
        prev_x = x;
        prev_z = z;
        prev_t = t;

        RungeKutta(dt, P, A);

        m_val = CalcMass(t);
        M_val = Mach(V);
        cR = Interpolate_cR(M_val);
        nR = nRalf(P, cR, M_val, A, m_val);
        alf = CalcAlpha(nza, nR);
        Cxa = Interpolate_Ca(alf, M_val);

        Predv << t << " "
            << m_val << " "
            << M_val << " "
            << alf << " "
            << cR << " "
            << Cxa << " "
            << V << " "
            << psi * 180 / PI << " "
            << x << " "
            << z << endl;

        if (!t_ayt_reached && t >= t_ayt) {
            double lambda = (t_ayt - prev_t) / dt;

            double V_ayt = prev_V + lambda * (V - prev_V);
            double psi_ayt = prev_psi + lambda * (psi - prev_psi);
            double x_ayt = prev_x + lambda * (x - prev_x);
            double z_ayt = prev_z + lambda * (z - prev_z);
            double M_ayt = V_ayt / a_atm;
            double m_ayt = CalcMass(t_ayt);
            double cR_ayt = Interpolate_cR(M_ayt);
            double nR_ayt = nRalf(P, cR_ayt, M_ayt, A, m_ayt);
            double alf_ayt = CalcAlpha(nza, nR_ayt);
            double Cxa_ayt = Interpolate_Ca(alf_ayt, M_ayt);

            Results << fixed << setprecision(4)
                << t_ayt << " "
                << m_ayt << " "
                << M_ayt << " "
                << alf_ayt << " "
                << cR_ayt << " "
                << Cxa_ayt << " "
                << V_ayt << " "
                << psi_ayt * 180 / PI << " "
                << x_ayt << " "
                << z_ayt << endl;

            cout << fixed << setprecision(2)
                << setw(8) << t_ayt
                << fixed << setprecision(4)
                << setw(10) << m_ayt
                << setw(10) << M_ayt
                << setw(12) << V_ayt
                << setw(10) << psi_ayt * 180 / PI
                << setw(12) << x_ayt
                << setw(12) << z_ayt
                << setw(10) << alf_ayt << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;

            t_ayt_reached = true;
        }

        if (abs(t - 0.5) < 0.001 || abs(t - 1.0) < 0.001 || abs(t - 1.5) < 0.001 ||
            abs(t - 2.0) < 0.001 || abs(t - 2.5) < 0.001 || abs(t - 3.0) < 0.001 ||
            abs(t - 3.5) < 0.001) {
            Results << fixed << setprecision(1)
                << t << " "
                << fixed << setprecision(4)
                << m_val << " "
                << M_val << " "
                << alf << " "
                << cR << " "
                << Cxa << " "
                << V << " "
                << psi * 180 / PI << " "
                << x << " "
                << z << endl;

            cout << fixed << setprecision(1)
                << setw(8) << t
                << fixed << setprecision(4)
                << setw(10) << m_val
                << setw(10) << M_val
                << setw(12) << V
                << setw(10) << psi * 180 / PI
                << setw(12) << x
                << setw(12) << z
                << setw(10) << alf << endl;
        }

        if (i == 43) {
            Results << fixed << setprecision(1)
                << t << " "
                << fixed << setprecision(4)
                << m_val << " "
                << M_val << " "
                << alf << " "
                << cR << " "
                << Cxa << " "
                << V << " "
                << psi * 180 / PI
                << " " << x
                << " " << z
                << endl;
        }
    }
    P = 0;

    bool t4_interpolated = false;
    int next_second = 5;

    while (M_val >= Mmin) {
        dt = 0.1;

        prev_V = V;
        prev_psi = psi;
        prev_x = x;
        prev_z = z;
        prev_t = t;

        RungeKutta(dt, P, A);

        m_val = CalcMass(t);
        M_val = Mach(V);
        cR = Interpolate_cR(M_val);
        nR = nRalf(P, cR, M_val, A, m_val);
        alf = CalcAlpha(nza, nR);
        Cxa = Interpolate_Ca(alf, M_val);

        Predv << fixed << setprecision(4)
            << t << " "
            << m_val << " "
            << M_val << " "
            << alf << " "
            << cR << " "
            << Cxa << " "
            << V << " "
            << psi * 180 / PI << " "
            << x << " "
            << z << endl;

        if (!t4_interpolated && prev_t < 4.0 && t >= 4.0) {
            double lambda = (4.0 - prev_t) / dt;

            double V_4 = prev_V + lambda * (V - prev_V);
            double psi_4 = prev_psi + lambda * (psi - prev_psi);
            double x_4 = prev_x + lambda * (x - prev_x);
            double z_4 = prev_z + lambda * (z - prev_z);
            double M_4 = V_4 / a_atm;
            double m_4 = CalcMass(4.0);
            double cR_4 = Interpolate_cR(M_4);
            double nR_4 = nRalf(P, cR_4, M_4, A, m_4);
            double alf_4 = CalcAlpha(nza, nR_4);
            double Cxa_4 = Interpolate_Ca(alf_4, M_4);

            Results << fixed << setprecision(1)
                << 4.0 << " "
                << fixed << setprecision(4)
                << m_4 << " "
                << M_4 << " "
                << alf_4 << " "
                << cR_4 << " "
                << Cxa_4 << " "
                << V_4 << " "
                << psi_4 * 180 / PI << " "
                << x_4 << " "
                << z_4 << endl;

            cout << fixed << setprecision(1)
                << setw(8) << 4.0
                << fixed << setprecision(4)
                << setw(10) << m_4
                << setw(10) << M_4
                << setw(12) << V_4
                << setw(10) << psi_4 * 180 / PI
                << setw(12) << x_4
                << setw(12) << z_4
                << setw(10) << alf_4 << endl;

            t4_interpolated = true;
        }

        if (abs(t - next_second) < 0.05) {
            Results << fixed << setprecision(1)
                << t << " "
                << fixed << setprecision(4)
                << m_val << " "
                << M_val << " "
                << alf << " "
                << cR << " "
                << Cxa << " "
                << V << " "
                << psi * 180 / PI << " "
                << x << " "
                << z << endl;

            cout << fixed << setprecision(1)
                << setw(8) << t
                << fixed << setprecision(4)
                << setw(10) << m_val
                << setw(10) << M_val
                << setw(12) << V
                << setw(10) << psi * 180 / PI
                << setw(12) << x
                << setw(12) << z
                << setw(10) << alf << endl;
            next_second++;
        }

        if (M_val < Mmin) {
            double lambda = (Vmin - prev_V) / (V - prev_V);
            double tmin = prev_t + lambda * dt;
            double psimin = prev_psi + lambda * (psi - prev_psi);
            double xmin = prev_x + lambda * (x - prev_x);
            double zmin = prev_z + lambda * (z - prev_z);

            double M_at_min = Mmin;
            double cR_at_min = Interpolate_cR(M_at_min);
            double nR_at_min = nRalf(P, cR_at_min, M_at_min, A, m_val);
            double alf_at_min = CalcAlpha(nza, nR_at_min);
            double Cxa_at_min = Interpolate_Ca(alf_at_min, M_at_min);

            Results << fixed << setprecision(4) << tmin << " " << m_val << " " << Mmin << " " << alf_at_min << " " << cR_at_min << " " << Cxa_at_min << " " << Vmin << " " << psimin * 180 / PI << " " << xmin << " " << zmin << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << fixed << setprecision(4) << setw(8) << tmin
                << setw(10) << m_val
                << setw(10) << Mmin
                << setw(12) << Vmin
                << setw(10) << psimin * 180 / PI
                << setw(12) << xmin
                << setw(12) << zmin
                << setw(10) << alf_at_min
                << endl;

            break;
        }
    }

    cout << "Результаты сохранены в папку проекта" << endl;
    return 0;
}