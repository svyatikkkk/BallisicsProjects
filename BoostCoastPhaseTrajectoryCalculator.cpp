#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

const double EarthRadius = 6356767.0;   // r        
const double GravitySea = 9.80665;      // g0
const double GasConstantAir = 287.053;  // R
const double AdiabaticIndex = 1.4;      // x
const double InitialMass = 1255;        // начальная масса 
const double MassFlux = 88;             // массовый расход
const double NozzleArea = 0.14;         // площадь сопла
const double PreasureSea = 101325;      // давление на уровне моря
const double MiddelArea = 0.231;        // площадь миделевого сечения 
const double MomentArm = 0.215;         // плечо момента
const double InertiaMoment = 127;       // момент инерции
const double InitialVelocity = 65;      // начальная скорость
const double InitialPitch = 0.820305;  // начальный угол тангажа              
const double ThrustSea = 2120;          // удельная тяга на уровне моря

// Структура для хранения данных
struct SimulationData {
    int N;
    double t;
    double m;
    double P;
    double V;
    double M;
    double Cxa;
    double alpha;
    double TetaC;
    double Cya;
    double dV_dt;
    double dTetaC_dt;
    double wz;
    double teta;
    double y;
    double dy_dt;
    double x;
    double dx_dt;
};

void CalcAtm(double h,
    double& H,
    double& T,
    double& P,
    double& p,
    double& a,
    double& g) {
    double T0 = 0, b = 0, H0 = 0, P0 = 0, Pl;
    g = GravitySea * pow(EarthRadius / (EarthRadius + h), 2);
    H = EarthRadius * h / (EarthRadius + h);

    if (H < 0) {
        H0 = -2000;
        T0 = 301.15;
        b = -0.0065;
        P0 = 127774;
    }

    if (H >= 0 && H < 11600) {
        H0 = 0;
        T0 = 288.150;
        b = -0.0065;
        P0 = 101325;
    }

    if (H >= 11600 && H < 20000) {
        H0 = 11600;
        T0 = 216.65;
        b = 0.0000;
        P0 = 20588.9;
    }

    if (H >= 20000 && H < 32000) {
        H0 = 20000;
        T0 = 216.65;
        b = 0.0010;
        P0 = 5474.87;
    }

    if (H >= 32000 && H < 47000) {
        H0 = 32000;
        T0 = 228.65;
        b = 0.0028;
        P0 = 868.014;
    }

    if (H >= 47000 && H < 51000) {
        H0 = 47000;
        T0 = 270.65;
        b = 0;
        P0 = 110.906;
    }

    if (H >= 51000 && H < 71000) {
        H0 = 51000;
        T0 = 270.65;
        b = -0.0028;
        P0 = 66.9384;
    }

    if (H >= 71000 && H < 85000) {
        H0 = 71000;
        T0 = 214.65;
        b = -0.002;
        P0 = 3.956390;
    }

    if (H >= 85000 && H < 94000) {
        H0 = 85000;
        T0 = 186.65;
        b = -0.002;
        P0 = 0.445710;
    }

    T = T0 + b * (H - H0);

    if (b != 0) {
        Pl = log10(P0) - GravitySea / (b * GasConstantAir) * log10(T / T0);
        P = pow(10, Pl);
    }
    else {
        Pl = log10(P0) - (0.434294 / (GasConstantAir * T)) * GravitySea * (H - H0);
        P = pow(10, Pl);
    }

    p = P / (GasConstantAir * T);

    a = sqrt(AdiabaticIndex * GasConstantAir * T);
}


void AlphaCalc(double& alpha, double tetaC, double teta) {
    alpha = teta - tetaC;
}


void MahCalc(double& Mah, double a_atm, double V) {
    Mah = V / a_atm;
}


void CxaCalc(double& Cxa, double Mah) {

    Cxa = 0.3;

    if (Mah >= 0.01 && Mah <= 0.55) {
        Cxa = 0.3;
    }

    if (Mah > 0.55 && Mah <= 1.1) {
        Cxa = 19.03 * pow(Mah, 3) - 7.877 * pow(Mah, 4) - 14.38 * pow(Mah, 2) + 4.125 * Mah - 0.0639;
    }

    if (Mah > 1.1 && Mah <= 10.0) {
        Cxa = 0.000572 * pow(Mah, 4) - 0.01499 * pow(Mah, 3) + 0.1415 * pow(Mah, 2) - 0.5836 * Mah + 1.358;
    }
}


void CyaCalc(double& Cya, double Mah) {

    Cya = 0.25;

    if (Mah >= 0.01 && Mah < 0.8) {
        Cya = 0.25;
    }

    if (Mah >= 0.8 && Mah < 0.9) {
        Cya = 0.65 - 0.5 * Mah;
    }

    if (Mah >= 0.9 && Mah < 1.0) {
        Cya = -0.7 + Mah;
    }

    if (Mah >= 1.0 && Mah < 1.06) {
        Cya = 0.14 + 0.16 * Mah;
    }

    if (Mah >= 1.06 && Mah < 1.1) {
        Cya = 1.9 - 1.5 * Mah;
    }

    if (Mah >= 1.1 && Mah <= 10.0) {
        Cya = 0.25;
    }
}


void ThrustCalc(double& thrust, double thrust_atm) {
    thrust = ThrustSea * MassFlux + NozzleArea * (PreasureSea - thrust_atm);
}


void XaCalc(double& Xa,
    double p_atm,
    double Cxa,
    double V) {

    double q;
    q = (p_atm / 2) * pow(V, 2);
    Xa = Cxa * MiddelArea * q;
}


void YaCalc(double& Ya,
    double p_atm,
    double Cya,
    double V,
    double alpha) {

    double q;
    q = (p_atm / 2) * pow(V, 2);
    Ya = Cya * alpha * MiddelArea * q;
}


void MassCalc(double& m, double t) {

    m = InitialMass - MassFlux * t;
}


void MzCalc(double& Mz,
    double Ya,
    double Xa,
    double alpha) {

    Mz = -(Ya * cos(alpha) + Xa * sin(alpha)) * MomentArm;
}


void GravForceXCalc(double& Fgx,
    double m,
    double g_atm,
    double tetaC) {

    Fgx = m * g_atm * sin(tetaC);
}


void GravForceYCalc(double& Fgy,
    double m,
    double g_atm,
    double TetaC) {

    Fgy = m * g_atm * cos(TetaC);
}


void dV_dtCalc(double& dV_dt,
    double P_atm,
    double p_atm,
    double Cxa,
    double V,
    double m,
    double g_atm,
    double tetaC,
    double t,
    double alpha) {

    double P, Xa, Fgx;
    double Fvar = 0;    //на будущее
    ThrustCalc(P, P_atm);
    XaCalc(Xa, p_atm, Cxa, V);
    GravForceXCalc(Fgx, m, g_atm, tetaC);

    dV_dt = (P * cos(alpha) - Xa - Fgx) / m + Fvar / m;
}


void dtetaC_dtCalc(double& dteta_dt,
    double P_atm,
    double p_atm,
    double Cya,
    double V,
    double alpha,
    double m,
    double g_atm,
    double tetaC,
    double t) {

    double P, Ya, Fgy;
    double Fvar = 0;    //на будущее
    ThrustCalc(P, P_atm);
    YaCalc(Ya, p_atm, Cya, V, alpha);
    GravForceYCalc(Fgy, m, g_atm, tetaC);

    dteta_dt = (P * sin(alpha) + Ya - Fgy) / (m * V) + Fvar / (m * V);
}


void dx_dtCalc(double& dx_dt,
    double V,
    double tetaC) {

    dx_dt = V * cos(tetaC);
}


void dy_dtCalc(double& dy_dt,
    double V,
    double tetaC) {

    dy_dt = V * sin(tetaC);
}


void dwz_dtCalc(double& dwz_dt, double Mz) {

    dwz_dt = Mz / InertiaMoment;
}


void dteta_dtCalc(double& dteta_dt, double wz) {

    dteta_dt = wz;
}


void EulerMethod_0(double& V_2, double& TetaC_2, double& x_2, double& y_2,
    double& wz_2, double& teta_2,
    double P_atm_1, double p_atm_1, double g_atm_1,
    double Cxa, double Cya, double V_1, double m,
    double TetaC_1, double t, double dt,
    double x_1, double y_1) {

    double dV_dt_1, dtetaC_dt_1, dx_dt, dy_dt;
    double alpha = 0.0;  // Для случая alpha = 0

    dV_dtCalc(dV_dt_1, P_atm_1, p_atm_1, Cxa, V_1, m, g_atm_1, TetaC_1, t, alpha);
    dtetaC_dtCalc(dtetaC_dt_1, P_atm_1, p_atm_1, Cya, V_1, alpha, m, g_atm_1, TetaC_1, t);
    dx_dtCalc(dx_dt, V_1, TetaC_1);
    dy_dtCalc(dy_dt, V_1, TetaC_1);

    V_2 = V_1 + dV_dt_1 * dt;
    TetaC_2 = TetaC_1 + dtetaC_dt_1 * dt;

    x_2 = x_1 + dx_dt * dt;
    y_2 = y_1 + dy_dt * dt;

    teta_2 = TetaC_2;  // При alpha = 0 углы совпадают
    wz_2 = 0;          // Угловая скорость = 0
}


void EulerMethod(double& V_2, double& TetaC_2, double& x_2, double& y_2,
    double& wz_2, double& teta_2, double P_atm_1, double p_atm_1,
    double g_atm_1, double Cxa, double Cya, double V_1, double m,
    double TetaC_1, double t, double alpha_1, double dt, double Mz,
    double x_1, double y_1, double wz_1, double teta_1) {

    double dV_dt_1, dtetaC_dt_1, dx_dt, dy_dt, dwz_dt, dteta_dt;

    dV_dtCalc(dV_dt_1, P_atm_1, p_atm_1, Cxa, V_1, m, g_atm_1, TetaC_1, t, alpha_1);
    dtetaC_dtCalc(dtetaC_dt_1, P_atm_1, p_atm_1, Cya, V_1, alpha_1, m, g_atm_1, TetaC_1, t);
    dx_dtCalc(dx_dt, V_1, TetaC_1);
    dy_dtCalc(dy_dt, V_1, TetaC_1);
    dwz_dtCalc(dwz_dt, Mz);
    dteta_dtCalc(dteta_dt, wz_1);

    V_2 = V_1 + dV_dt_1 * dt;
    TetaC_2 = TetaC_1 + dtetaC_dt_1 * dt;

    x_2 = x_1 + dx_dt * dt;
    y_2 = y_1 + dy_dt * dt;

    teta_2 = teta_1 + dteta_dt * dt;
    wz_2 = wz_1 + dwz_dt * dt;
}

// Модифицированный метод Эйлера для alpha = 0
void EulerModMethod_0(double& dV_dt_1, double& dTetaC_dt_1, double& dx_dt_1, double& dy_dt_1,
    double& V_2, double& TetaC_2, double& x_2, double& y_2,
    double& teta_2, double& wz_2, double dt,
    double x_1, double y_1, double Cxa, double V_1,
    double m, double TetaC_1, double t, double Cya) {

    double dV, dTetaC, dx, dy, dV_dt, dTetaC_dt, dx_dt, dy_dt;
    double alpha = 0.0;  // Для случая alpha = 0
    double Mah_1, Cxa_1, Cya_1;
    double H_atm, T_atm, P_atm, p_atm, a_atm, g_atm;
    double H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1;

    CalcAtm(y_1, H_atm, T_atm, P_atm, p_atm, a_atm, g_atm);

    dV_dtCalc(dV_dt, P_atm, p_atm, Cxa, V_1, m, g_atm, TetaC_1, t, alpha);
    dtetaC_dtCalc(dTetaC_dt, P_atm, p_atm, Cya, V_1, alpha, m, g_atm, TetaC_1, t);
    dx_dtCalc(dx_dt, V_1, TetaC_1);
    dy_dtCalc(dy_dt, V_1, TetaC_1);

    dV = dV_dt * dt;
    dTetaC = dTetaC_dt * dt;
    dx = dx_dt * dt;
    dy = dy_dt * dt;

    CalcAtm(y_1 + dy / 2, H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1);
    MahCalc(Mah_1, a_atm_1, V_1 + dV / 2);
    CxaCalc(Cxa_1, Mah_1);
    CyaCalc(Cya_1, Mah_1);

    dV_dtCalc(dV_dt_1, P_atm_1, p_atm_1, Cxa_1, V_1 + dV / 2, m, g_atm_1, TetaC_1 + dTetaC / 2, t + dt / 2, alpha);
    dtetaC_dtCalc(dTetaC_dt_1, P_atm_1, p_atm_1, Cya_1, V_1 + dV / 2, alpha, m, g_atm_1, TetaC_1 + dTetaC / 2, t + dt / 2);
    dx_dtCalc(dx_dt_1, V_1 + dV / 2, TetaC_1 + dTetaC / 2);
    dy_dtCalc(dy_dt_1, V_1 + dV / 2, TetaC_1 + dTetaC / 2);

    V_2 = V_1 + dV_dt_1 * dt;
    TetaC_2 = TetaC_1 + dTetaC_dt_1 * dt;

    x_2 = x_1 + dx_dt_1 * dt;
    y_2 = y_1 + dy_dt_1 * dt;

    teta_2 = TetaC_2;  // При alpha = 0 углы совпадают
    wz_2 = 0;          // Угловая скорость = 0
}

// Модифицированный метод Эйлера для alpha != 0
void EulerModMethod(double& dV_dt_1, double& dwz_dt_1, double& dteta_dt_1,
    double& dTetaC_dt_1, double& dx_dt_1, double& dy_dt_1,
    double& V_2, double& TetaC_2, double& x_2, double& y_2,
    double& teta_2, double& wz_2, double dt, double wz_1,
    double Mz, double x_1, double y_1, double Cxa, double V_1,
    double m, double TetaC_1, double t, double alpha,
    double Cya, double teta_1) {

    double dV, dTetaC, dx, dy, dteta, dwz, dV_dt, dTetaC_dt, dx_dt, dy_dt, dteta_dt, dwz_dt, alpha_1, Mah_1, Cxa_1, Cya_1, Mz_1, Xa_1, Ya_1, H_atm, T_atm, P_atm, p_atm, a_atm, g_atm, H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1;

    CalcAtm(y_1, H_atm, T_atm, P_atm, p_atm, a_atm, g_atm);

    dV_dtCalc(dV_dt, P_atm, p_atm, Cxa, V_1, m, g_atm, TetaC_1, t, alpha);
    dtetaC_dtCalc(dTetaC_dt, P_atm, p_atm, Cya, V_1, alpha, m, g_atm, TetaC_1, t);
    dx_dtCalc(dx_dt, V_1, TetaC_1);
    dy_dtCalc(dy_dt, V_1, TetaC_1);
    dteta_dtCalc(dteta_dt, wz_1);
    dwz_dtCalc(dwz_dt, Mz);

    dV = dV_dt * dt;
    dTetaC = dTetaC_dt * dt;
    dx = dx_dt * dt;
    dy = dy_dt * dt;
    dteta = dteta_dt * dt;;
    dwz = dwz_dt * dt;

    CalcAtm(y_1 + dy / 2, H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1);
    AlphaCalc(alpha_1, TetaC_1 + dTetaC / 2, teta_1 + dteta / 2);
    MahCalc(Mah_1, a_atm_1, V_1 + dV / 2);
    CxaCalc(Cxa_1, Mah_1);
    CyaCalc(Cya_1, Mah_1);
    XaCalc(Xa_1, p_atm_1, Cxa_1, V_1 + dV / 2);
    YaCalc(Ya_1, p_atm_1, Cya_1, V_1 + dV / 2, alpha_1);
    MzCalc(Mz_1, Ya_1, Xa_1, alpha_1);

    dV_dtCalc(dV_dt_1, P_atm_1, p_atm_1, Cxa_1, V_1 + dV / 2, m, g_atm_1, TetaC_1 + dTetaC / 2, t + dt / 2, alpha_1);
    dtetaC_dtCalc(dTetaC_dt_1, P_atm_1, p_atm_1, Cya_1, V_1 + dV / 2, alpha_1, m, g_atm_1, TetaC_1 + dTetaC / 2, t + dt / 2);
    dx_dtCalc(dx_dt_1, V_1 + dV / 2, TetaC_1 + dTetaC / 2);
    dy_dtCalc(dy_dt_1, V_1 + dV / 2, TetaC_1 + dTetaC / 2);
    dteta_dtCalc(dteta_dt_1, wz_1 + dwz / 2);
    dwz_dtCalc(dwz_dt_1, Mz_1);

    V_2 = V_1 + dV_dt_1 * dt;
    TetaC_2 = TetaC_1 + dTetaC_dt_1 * dt;

    x_2 = x_1 + dx_dt_1 * dt;
    y_2 = y_1 + dy_dt_1 * dt;

    teta_2 = teta_1 + dteta_dt_1 * dt;
    wz_2 = wz_1 + dwz_dt_1 * dt;
}

// Метод Рунге-Кутта для alpha = 0
void RungeKuttMethod_0(double& V_2, double& TetaC_2, double& x_2, double& y_2,
    double& wz_2, double& teta_2, double dt, double y_1,
    double V_1, double m, double TetaC_1, double t,
    double Cxa, double Cya, double x_1) {

    double H_atm, T_atm, P_atm, p_atm, a_atm, g_atm, dV_dt, dTetaC_dt, dx_dt, dy_dt;
    double k1_V, k1_TetaC, k1_x, k1_y;
    double H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1;
    double alpha = 0.0, Mah_1, Cxa_1, Cya_1, dV_dt_1, dTetaC_dt_1, dx_dt_1, dy_dt_1;
    double k2_V, k2_TetaC, k2_x, k2_y, H_atm_2, T_atm_2, P_atm_2, p_atm_2, a_atm_2, g_atm_2;
    double Mah_2, Cxa_2, Cya_2, dV_dt_2, dTetaC_dt_2, dx_dt_2, dy_dt_2;
    double k3_V, k3_TetaC, k3_x, k3_y, H_atm_3, T_atm_3, P_atm_3, p_atm_3, a_atm_3, g_atm_3;
    double Mah_3, Cxa_3, Cya_3, dV_dt_3, dTetaC_dt_3, dx_dt_3, dy_dt_3;
    double k4_V, k4_TetaC, k4_x, k4_y;

    CalcAtm(y_1, H_atm, T_atm, P_atm, p_atm, a_atm, g_atm);

    dV_dtCalc(dV_dt, P_atm, p_atm, Cxa, V_1, m, g_atm, TetaC_1, t, alpha);
    dtetaC_dtCalc(dTetaC_dt, P_atm, p_atm, Cya, V_1, alpha, m, g_atm, TetaC_1, t);
    dx_dtCalc(dx_dt, V_1, TetaC_1);
    dy_dtCalc(dy_dt, V_1, TetaC_1);

    k1_V = dt * dV_dt;
    k1_TetaC = dt * dTetaC_dt;
    k1_x = dt * dx_dt;
    k1_y = dt * dy_dt;

    CalcAtm(y_1 + k1_y / 2, H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1);
    MahCalc(Mah_1, a_atm_1, V_1 + k1_V / 2);
    CxaCalc(Cxa_1, Mah_1);
    CyaCalc(Cya_1, Mah_1);

    dV_dtCalc(dV_dt_1, P_atm_1, p_atm_1, Cxa_1, V_1 + k1_V / 2, m, g_atm_1, TetaC_1 + k1_TetaC / 2, t + dt / 2, alpha);
    dtetaC_dtCalc(dTetaC_dt_1, P_atm_1, p_atm_1, Cya_1, V_1 + k1_V / 2, alpha, m, g_atm_1, TetaC_1 + k1_TetaC / 2, t + dt / 2);
    dx_dtCalc(dx_dt_1, V_1 + k1_V / 2, TetaC_1 + k1_TetaC / 2);
    dy_dtCalc(dy_dt_1, V_1 + k1_V / 2, TetaC_1 + k1_TetaC / 2);

    k2_V = dt * dV_dt_1;
    k2_TetaC = dt * dTetaC_dt_1;
    k2_x = dt * dx_dt_1;
    k2_y = dt * dy_dt_1;

    CalcAtm(y_1 + k2_y / 2, H_atm_2, T_atm_2, P_atm_2, p_atm_2, a_atm_2, g_atm_2);
    MahCalc(Mah_2, a_atm_2, V_1 + k2_V / 2);
    CxaCalc(Cxa_2, Mah_2);
    CyaCalc(Cya_2, Mah_2);

    dV_dtCalc(dV_dt_2, P_atm_2, p_atm_2, Cxa_2, V_1 + k2_V / 2, m, g_atm_2, TetaC_1 + k2_TetaC / 2, t + dt / 2, alpha);
    dtetaC_dtCalc(dTetaC_dt_2, P_atm_2, p_atm_2, Cya_2, V_1 + k2_V / 2, alpha, m, g_atm_2, TetaC_1 + k2_TetaC / 2, t + dt / 2);
    dx_dtCalc(dx_dt_2, V_1 + k2_V / 2, TetaC_1 + k2_TetaC / 2);
    dy_dtCalc(dy_dt_2, V_1 + k2_V / 2, TetaC_1 + k2_TetaC / 2);

    k3_V = dt * dV_dt_2;
    k3_TetaC = dt * dTetaC_dt_2;
    k3_x = dt * dx_dt_2;
    k3_y = dt * dy_dt_2;

    CalcAtm(y_1 + k3_y, H_atm_3, T_atm_3, P_atm_3, p_atm_3, a_atm_3, g_atm_3);
    MahCalc(Mah_3, a_atm_3, V_1 + k3_V);
    CxaCalc(Cxa_3, Mah_3);
    CyaCalc(Cya_3, Mah_3);

    dV_dtCalc(dV_dt_3, P_atm_3, p_atm_3, Cxa_3, V_1 + k3_V, m, g_atm_3, TetaC_1 + k3_TetaC, t + dt, alpha);
    dtetaC_dtCalc(dTetaC_dt_3, P_atm_3, p_atm_3, Cya_3, V_1 + k3_V, alpha, m, g_atm_3, TetaC_1 + k3_TetaC, t + dt);
    dx_dtCalc(dx_dt_3, V_1 + k3_V, TetaC_1 + k3_TetaC);
    dy_dtCalc(dy_dt_3, V_1 + k3_V, TetaC_1 + k3_TetaC);

    k4_V = dt * dV_dt_3;
    k4_TetaC = dt * dTetaC_dt_3;
    k4_x = dt * dx_dt_3;
    k4_y = dt * dy_dt_3;

    V_2 = V_1 + (k1_V + 2 * k2_V + 2 * k3_V + k4_V) / 6;
    TetaC_2 = TetaC_1 + (k1_TetaC + 2 * k2_TetaC + 2 * k3_TetaC + k4_TetaC) / 6;
    x_2 = x_1 + (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    y_2 = y_1 + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;

    teta_2 = TetaC_2;  // При alpha = 0 углы совпадают
    wz_2 = 0;          // Угловая скорость = 0
}

// Метод Рунге-Кутта для alpha != 0
void RungeKuttMethod(double& V_2, double& TetaC_2, double& x_2, double& y_2,
    double& wz_2, double& teta_2, double dt, double y_1,
    double V_1, double m, double TetaC_1, double t,
    double teta_1, double wz_1, double Mz, double Cxa,
    double Cya, double alpha_1, double x_1) {

    double H_atm, T_atm, P_atm, p_atm, a_atm, g_atm, dV_dt, dTetaC_dt, dx_dt, dy_dt, dteta_dt, dwz_dt;
    double k1_V, k1_TetaC, k1_x, k1_y, k1_teta, k1_wz;
    double H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1;
    double alpha_2, Mah_1, Cxa_1, Cya_1, Xa_1, Ya_1, Mz_1, dV_dt_1, dTetaC_dt_1, dx_dt_1, dy_dt_1, dteta_dt_1, dwz_dt_1;
    double k2_V, k2_TetaC, k2_x, k2_y, k2_teta, k2_wz, H_atm_2, T_atm_2, P_atm_2, p_atm_2, a_atm_2, g_atm_2;
    double alpha_3, Mah_2, Cxa_2, Cya_2, Xa_2, Ya_2, Mz_2, dV_dt_2, dTetaC_dt_2, dx_dt_2, dy_dt_2, dteta_dt_2, dwz_dt_2;
    double k3_V, k3_TetaC, k3_x, k3_y, k3_teta, k3_wz, H_atm_3, T_atm_3, P_atm_3, p_atm_3, a_atm_3, g_atm_3;
    double alpha_4, Mah_3, Cxa_3, Cya_3, Xa_3, Ya_3, Mz_3, dV_dt_3, dTetaC_dt_3, dx_dt_3, dy_dt_3, dteta_dt_3, dwz_dt_3;
    double k4_V, k4_TetaC, k4_x, k4_y, k4_teta, k4_wz;

    CalcAtm(y_1, H_atm, T_atm, P_atm, p_atm, a_atm, g_atm);

    dV_dtCalc(dV_dt, P_atm, p_atm, Cxa, V_1, m, g_atm, TetaC_1, t, alpha_1);
    dtetaC_dtCalc(dTetaC_dt, P_atm, p_atm, Cya, V_1, alpha_1, m, g_atm, TetaC_1, t);
    dx_dtCalc(dx_dt, V_1, TetaC_1);
    dy_dtCalc(dy_dt, V_1, TetaC_1);
    dteta_dtCalc(dteta_dt, wz_1);
    dwz_dtCalc(dwz_dt, Mz);

    k1_V = dt * dV_dt;
    k1_TetaC = dt * dTetaC_dt;
    k1_x = dt * dx_dt;
    k1_y = dt * dy_dt;
    k1_teta = dt * dteta_dt;
    k1_wz = dt * dwz_dt;

    CalcAtm(y_1 + k1_y / 2, H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1);
    AlphaCalc(alpha_2, TetaC_1 + k1_TetaC / 2, teta_1 + k1_teta / 2);
    MahCalc(Mah_1, a_atm_1, V_1 + k1_V / 2);
    CxaCalc(Cxa_1, Mah_1);
    CyaCalc(Cya_1, Mah_1);
    XaCalc(Xa_1, p_atm_1, Cxa_1, V_1 + k1_V / 2);
    YaCalc(Ya_1, p_atm_1, Cya_1, V_1 + k1_V / 2, alpha_2);
    MzCalc(Mz_1, Ya_1, Xa_1, alpha_2);

    dV_dtCalc(dV_dt_1, P_atm_1, p_atm_1, Cxa_1, V_1 + k1_V / 2, m, g_atm_1, TetaC_1 + k1_TetaC / 2, t + dt / 2, alpha_2);
    dtetaC_dtCalc(dTetaC_dt_1, P_atm_1, p_atm_1, Cya_1, V_1 + k1_V / 2, alpha_2, m, g_atm_1, TetaC_1 + k1_TetaC / 2, t + dt / 2);
    dx_dtCalc(dx_dt_1, V_1 + k1_V / 2, TetaC_1 + k1_TetaC / 2);
    dy_dtCalc(dy_dt_1, V_1 + k1_V / 2, TetaC_1 + k1_TetaC / 2);
    dteta_dtCalc(dteta_dt_1, wz_1 + k1_wz / 2);
    dwz_dtCalc(dwz_dt_1, Mz_1);

    k2_V = dt * dV_dt_1;
    k2_TetaC = dt * dTetaC_dt_1;
    k2_x = dt * dx_dt_1;
    k2_y = dt * dy_dt_1;
    k2_teta = dt * dteta_dt_1;
    k2_wz = dt * dwz_dt_1;

    CalcAtm(y_1 + k2_y / 2, H_atm_2, T_atm_2, P_atm_2, p_atm_2, a_atm_2, g_atm_2);
    AlphaCalc(alpha_3, TetaC_1 + k2_TetaC / 2, teta_1 + k2_teta / 2);
    MahCalc(Mah_2, a_atm_2, V_1 + k2_V / 2);
    CxaCalc(Cxa_2, Mah_2);
    CyaCalc(Cya_2, Mah_2);
    XaCalc(Xa_2, p_atm_2, Cxa_2, V_1 + k2_V / 2);
    YaCalc(Ya_2, p_atm_2, Cya_2, V_1 + k2_V / 2, alpha_3);
    MzCalc(Mz_2, Ya_2, Xa_2, alpha_3);

    dV_dtCalc(dV_dt_2, P_atm_2, p_atm_2, Cxa_2, V_1 + k2_V / 2, m, g_atm_2, TetaC_1 + k2_TetaC / 2, t + dt / 2, alpha_3);
    dtetaC_dtCalc(dTetaC_dt_2, P_atm_2, p_atm_2, Cya_2, V_1 + k2_V / 2, alpha_3, m, g_atm_2, TetaC_1 + k2_TetaC / 2, t + dt / 2);
    dx_dtCalc(dx_dt_2, V_1 + k2_V / 2, TetaC_1 + k2_TetaC / 2);
    dy_dtCalc(dy_dt_2, V_1 + k2_V / 2, TetaC_1 + k2_TetaC / 2);
    dteta_dtCalc(dteta_dt_2, wz_1 + k2_wz / 2);
    dwz_dtCalc(dwz_dt_2, Mz_2);

    k3_V = dt * dV_dt_2;
    k3_TetaC = dt * dTetaC_dt_2;
    k3_x = dt * dx_dt_2;
    k3_y = dt * dy_dt_2;
    k3_teta = dt * dteta_dt_2;
    k3_wz = dt * dwz_dt_2;

    CalcAtm(y_1 + k3_y, H_atm_3, T_atm_3, P_atm_3, p_atm_3, a_atm_3, g_atm_3);
    AlphaCalc(alpha_4, TetaC_1 + k3_TetaC, teta_1 + k3_teta);
    MahCalc(Mah_3, a_atm_3, V_1 + k3_V);
    CxaCalc(Cxa_3, Mah_3);
    CyaCalc(Cya_3, Mah_3);
    XaCalc(Xa_3, p_atm_3, Cxa_3, V_1 + k3_V);
    YaCalc(Ya_3, p_atm_3, Cya_3, V_1 + k3_V, alpha_4);
    MzCalc(Mz_3, Ya_3, Xa_3, alpha_4);

    dV_dtCalc(dV_dt_3, P_atm_3, p_atm_3, Cxa_3, V_1 + k3_V, m, g_atm_3, TetaC_1 + k3_TetaC, t + dt, alpha_4);
    dtetaC_dtCalc(dTetaC_dt_3, P_atm_3, p_atm_3, Cya_3, V_1 + k3_V, alpha_4, m, g_atm_3, TetaC_1 + k3_TetaC, t + dt);
    dx_dtCalc(dx_dt_3, V_1 + k3_V, TetaC_1 + k3_TetaC);
    dy_dtCalc(dy_dt_3, V_1 + k3_V, TetaC_1 + k3_TetaC);
    dteta_dtCalc(dteta_dt_3, wz_1 + k3_wz);
    dwz_dtCalc(dwz_dt_3, Mz_3);

    k4_V = dt * dV_dt_3;
    k4_TetaC = dt * dTetaC_dt_3;
    k4_x = dt * dx_dt_3;
    k4_y = dt * dy_dt_3;
    k4_teta = dt * dteta_dt_3;
    k4_wz = dt * dwz_dt_3;

    V_2 = V_1 + (k1_V + 2 * k2_V + 2 * k3_V + k4_V) / 6;
    TetaC_2 = TetaC_1 + (k1_TetaC + 2 * k2_TetaC + 2 * k3_TetaC + k4_TetaC) / 6;
    x_2 = x_1 + (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    y_2 = y_1 + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;

    teta_2 = teta_1 + (k1_teta + 2 * k2_teta + 2 * k3_teta + k4_teta) / 6;
    wz_2 = wz_1 + (k1_wz + 2 * k2_wz + 2 * k3_wz + k4_wz) / 6;
}

// Функция для записи данных в CSV файл
void writeToCSV(const string& filename, const vector<SimulationData>& data) {
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Записываем заголовки
    file << "N;t;m;P;V;M;Cxa;alpha;TetaC;Cya;dV_dt;dTetaC_dt;wz;teta;y;dy_dt;x;dx_dt" << endl;

    // Записываем данные с разделителем ";" (точка с запятой)
    for (const auto& row : data) {
        file << fixed << setprecision(6);
        file << row.N << ";"
            << row.t << ";"
            << row.m << ";"
            << row.P << ";"
            << row.V << ";"
            << row.M << ";"
            << row.Cxa << ";"
            << row.alpha << ";"
            << row.TetaC << ";"
            << row.Cya << ";"
            << row.dV_dt << ";"
            << row.dTetaC_dt << ";"
            << row.wz << ";"
            << row.teta << ";"
            << row.y << ";"
            << row.dy_dt << ";"
            << row.x << ";"
            << row.dx_dt << endl;
    }

    file.close();
    cout << "Data written to " << filename << endl;
}

// Общая функция для выполнения симуляции
void runSimulation(const string& methodName, int K, double dt, int methodType) {
    double t = 0, InitialVelocity = 65,
        TetaC_0 = 0.820305, teta_0 = 0.820305,
        wz_0 = 0.03, x0 = 0, y0 = 2825;
    double V_1, TetaC_1, teta_1, x_1, y_1, wz_1, y_2, TetaC_2,
        teta_2, V_2, x_2, wz_2;
    double Mah_1, Cxa_1, Cya_1, alpha_1, P, Ya, Xa, Mz,
        m, dV_dt, dTetaC_dt_1, dx_dt, dy_dt;
    double H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1;
    int N = 0;
    int iteration = 0;  // Счетчик всех итераций

    // Определяем шаг вывода в таблицу
    int output_step = 1;
    if (fabs(dt - 0.01) < 1e-10) {
        output_step = 10;  // Для dt=0.01 выводим каждую 10-ю строку
    }
    else if (fabs(dt - 0.001) < 1e-10) {
        output_step = 100; // Для dt=0.001 выводим каждую 100-ю строку
    }

    vector<SimulationData> allData;
    // Вектор для хранения последних двух точек для интерполяции
    vector<SimulationData> lastTwoPoints;

    // Инициализация переменных
    V_1 = InitialVelocity;
    V_2 = InitialVelocity;
    TetaC_1 = TetaC_0;
    TetaC_2 = TetaC_0;
    teta_1 = teta_0;
    teta_2 = teta_0;
    x_1 = x0;
    x_2 = x0;
    y_1 = y0;
    y_2 = y0;
    wz_1 = wz_0;
    wz_2 = wz_0;
    m = InitialMass;

    cout << "\n=== " << methodName << " ===" << endl;
    cout << "dt = " << dt << ", K = " << K << endl;
    cout << "Output step: every " << output_step << " iterations" << endl;

    // Основной цикл до 2.4 включительно
    while (t <= 2.4 + 1e-10) {
        // Расчет атмосферы
        CalcAtm(y_1, H_atm_1, T_atm_1, P_atm_1, p_atm_1, a_atm_1, g_atm_1);

        // Расчет угла атаки
        if (K == 0) {
            alpha_1 = 0.0;
        }
        else {
            AlphaCalc(alpha_1, TetaC_1, teta_1);
        }

        MahCalc(Mah_1, a_atm_1, V_1);
        CxaCalc(Cxa_1, Mah_1);
        CyaCalc(Cya_1, Mah_1);
        ThrustCalc(P, P_atm_1);

        if (t > 2.49) {
            P = 0;
        }

        if (K == 0) {
            Ya = 0;
        }
        else {
            YaCalc(Ya, p_atm_1, Cya_1, V_1, alpha_1);
        }

        XaCalc(Xa, p_atm_1, Cxa_1, V_1);

        if (K == 0) {
            Mz = 0;
        }
        else {
            MzCalc(Mz, Ya, Xa, alpha_1);
        }

        // Расчет массы
        if (t <= 2.49) {
            MassCalc(m, t);
        }
        else {
            m = InitialMass - MassFlux * 2.49;
        }

        dV_dtCalc(dV_dt, P_atm_1, p_atm_1, Cxa_1, V_1, m, g_atm_1, TetaC_1, t, alpha_1);
        dtetaC_dtCalc(dTetaC_dt_1, P_atm_1, p_atm_1, Cya_1, V_1, alpha_1, m, g_atm_1, TetaC_1, t);
        dx_dtCalc(dx_dt, V_1, TetaC_1);
        dy_dtCalc(dy_dt, V_1, TetaC_1);

        // Выбор метода интегрирования
        if (methodType == 1) { // Метод Эйлера
            if (K == 0) {
                EulerMethod_0(V_2, TetaC_2, x_2, y_2, wz_2, teta_2,
                    P_atm_1, p_atm_1, g_atm_1,
                    Cxa_1, Cya_1, V_1, m, TetaC_1, t, dt,
                    x_1, y_1);
            }
            else {
                EulerMethod(V_2, TetaC_2, x_2, y_2, wz_2, teta_2,
                    P_atm_1, p_atm_1, g_atm_1,
                    Cxa_1, Cya_1, V_1, m, TetaC_1, t, alpha_1, dt, Mz,
                    x_1, y_1, wz_1, teta_1);
            }
        }
        else if (methodType == 2) { // Модифицированный метод Эйлера
            double dV_dt_1, dx_dt_1, dy_dt_1, dwz_dt_1 = 0, dteta_dt_1 = 0, dTetaC_dt_1_var;

            if (K == 0) {
                EulerModMethod_0(dV_dt_1, dTetaC_dt_1_var, dx_dt_1, dy_dt_1,
                    V_2, TetaC_2, x_2, y_2, teta_2, wz_2,
                    dt, x_1, y_1, Cxa_1, V_1, m, TetaC_1, t, Cya_1);
            }
            else {
                EulerModMethod(dV_dt_1, dwz_dt_1, dteta_dt_1, dTetaC_dt_1_var,
                    dx_dt_1, dy_dt_1, V_2, TetaC_2, x_2, y_2,
                    teta_2, wz_2, dt, wz_1, Mz, x_1, y_1,
                    Cxa_1, V_1, m, TetaC_1, t, alpha_1, Cya_1, teta_1);
            }
        }
        else if (methodType == 3) { // Метод Рунге-Кутта
            if (K == 0) {
                RungeKuttMethod_0(V_2, TetaC_2, x_2, y_2, wz_2, teta_2,
                    dt, y_1, V_1, m, TetaC_1, t,
                    Cxa_1, Cya_1, x_1);
            }
            else {
                RungeKuttMethod(V_2, TetaC_2, x_2, y_2, wz_2, teta_2,
                    dt, y_1, V_1, m, TetaC_1, t, teta_1, wz_1,
                    Mz, Cxa_1, Cya_1, alpha_1, x_1);
            }
        }

        // Сохраняем данные только для каждой output_step-ой итерации
        if (iteration % output_step == 0) {
            SimulationData data;
            data.N = N;  // Номер записи в таблице
            data.t = t;
            data.m = m;
            data.P = P;
            data.V = V_1;
            data.M = Mah_1;
            data.Cxa = Cxa_1;
            data.alpha = alpha_1 * 180 / 3.141592653;  // Преобразуем в градусы
            data.TetaC = TetaC_1 * 180 / 3.141592653;  // Преобразуем в градусы
            data.Cya = Cya_1;
            data.dV_dt = dV_dt;
            data.dTetaC_dt = dTetaC_dt_1;
            data.wz = wz_1;
            data.teta = teta_1 * 180 / 3.141592653;  // Преобразуем в градусы
            data.y = y_1;
            data.dy_dt = dy_dt;
            data.x = x_1;
            data.dx_dt = dx_dt;

            allData.push_back(data);

            // Сохраняем последние две точки для интерполяции
            lastTwoPoints.push_back(data);
            if (lastTwoPoints.size() > 2) {
                lastTwoPoints.erase(lastTwoPoints.begin());
            }

            N++;  // Увеличиваем счетчик записей в таблице
        }

        // Вывод в консоль каждые 100*output_step итераций
        if (iteration % (100 * output_step) == 0) {
            cout << "Iteration: " << iteration << ", t = " << t << ", V = " << V_1 << ", y = " << y_1 << endl;
        }

        t += dt;
        iteration++;

        // Обновление переменных для следующей итерации
        V_1 = V_2;
        TetaC_1 = TetaC_2;
        teta_1 = teta_2;
        x_1 = x_2;
        y_1 = y_2;
        wz_1 = wz_2;

        // Если достигли 2.4, выходим из цикла
        if (t > 2.4 + 1e-10) {
            break;
        }
    }

    // Теперь добавляем точку 2.49 с помощью линейной интерполяции
    if (lastTwoPoints.size() >= 2) {
        SimulationData point_2_3 = lastTwoPoints[0];  // Точка при t=2.3
        SimulationData point_2_4 = lastTwoPoints[1];  // Точка при t=2.4

        // Коэффициент интерполяции: (2.49 - 2.3) / (2.4 - 2.3) = 1.9
        double t1 = 2.3;
        double t2 = 2.4;
        double t_target = 2.49;
        double alpha = (t_target - t1) / (t2 - t1);  // = 1.9

        SimulationData interpolated_data;
        interpolated_data.N = N;
        interpolated_data.t = t_target;

        // Линейная интерполяция всех величин
        interpolated_data.m = point_2_3.m + alpha * (point_2_4.m - point_2_3.m);
        interpolated_data.P = point_2_3.P + alpha * (point_2_4.P - point_2_3.P);
        interpolated_data.V = point_2_3.V + alpha * (point_2_4.V - point_2_3.V);
        interpolated_data.M = point_2_3.M + alpha * (point_2_4.M - point_2_3.M);
        interpolated_data.Cxa = point_2_3.Cxa + alpha * (point_2_4.Cxa - point_2_3.Cxa);
        interpolated_data.alpha = point_2_3.alpha + alpha * (point_2_4.alpha - point_2_3.alpha);
        interpolated_data.TetaC = point_2_3.TetaC + alpha * (point_2_4.TetaC - point_2_3.TetaC);
        interpolated_data.Cya = point_2_3.Cya + alpha * (point_2_4.Cya - point_2_3.Cya);
        interpolated_data.dV_dt = point_2_3.dV_dt + alpha * (point_2_4.dV_dt - point_2_3.dV_dt);
        interpolated_data.dTetaC_dt = point_2_3.dTetaC_dt + alpha * (point_2_4.dTetaC_dt - point_2_3.dTetaC_dt);
        interpolated_data.wz = point_2_3.wz + alpha * (point_2_4.wz - point_2_3.wz);
        interpolated_data.teta = point_2_3.teta + alpha * (point_2_4.teta - point_2_3.teta);
        interpolated_data.y = point_2_3.y + alpha * (point_2_4.y - point_2_3.y);
        interpolated_data.dy_dt = point_2_3.dy_dt + alpha * (point_2_4.dy_dt - point_2_3.dy_dt);
        interpolated_data.x = point_2_3.x + alpha * (point_2_4.x - point_2_3.x);
        interpolated_data.dx_dt = point_2_3.dx_dt + alpha * (point_2_4.dx_dt - point_2_3.dx_dt);

        allData.push_back(interpolated_data);

        cout << "\nAdded interpolated point at t=2.49:" << endl;
        cout << "Interpolated from t1=" << t1 << " to t2=" << t2 << " with alpha=" << alpha << endl;
        cout << "V = " << interpolated_data.V << ", y = " << interpolated_data.y << ", x = " << interpolated_data.x << endl;
    }

    // Создаем имя файла
    string filename = methodName;
    if (K == 0) {
        filename += "_alpha0";
    }
    else {
        filename += "_alpha";
    }
    // Добавляем информацию о шаге в имя файла
    if (fabs(dt - 0.1) < 1e-10) {
        filename += "_dt01";
    }
    else if (fabs(dt - 0.01) < 1e-10) {
        filename += "_dt001";
    }
    else if (fabs(dt - 0.001) < 1e-10) {
        filename += "_dt0001";
    }
    else {
        filename += "_dt" + to_string(dt).substr(0, 4);
    }
    filename += ".csv";

    // Записываем данные в CSV файл
    writeToCSV(filename, allData);

    cout << "\nSimulation completed." << endl;
    cout << "Total iterations: " << iteration << endl;
    cout << "Data points in table: " << allData.size() << endl;
    cout << "Data saved to: " << filename << endl;
}
int main() {
    int metod, K;
    double dt;

    cout << "Choose numerical method:" << endl;
    cout << "1 - Euler's method" << endl;
    cout << "2 - Modified Euler's method" << endl;
    cout << "3 - Runge-Kutta method" << endl;
    cout << "4 - Run ALL methods (comparison)" << endl;
    cin >> metod;

    cout << "Choose angle of attack:" << endl;
    cout << "0 - alpha = 0" << endl;
    cout << "1 - alpha != 0" << endl;
    cin >> K;

    cout << "Input time step dt: ";
    cin >> dt;

    if (metod == 1) {
        runSimulation("Euler_Method", K, dt, 1);
    }
    else if (metod == 2) {
        runSimulation("Modified_Euler_Method", K, dt, 2);
    }
    else if (metod == 3) {
        runSimulation("Runge_Kutta_Method", K, dt, 3);
    }
    else if (metod == 4) {
        // Запускаем все методы для сравнения
        cout << "\n=== Running ALL methods for comparison ===" << endl;
        runSimulation("Euler_Method", K, dt, 1);
        runSimulation("Modified_Euler_Method", K, dt, 2);
        runSimulation("Runge_Kutta_Method", K, dt, 3);
        cout << "\nAll simulations completed! Data saved to files." << endl;
    }
    else {
        cout << "Invalid method choice!" << endl;
        return 1;
    }

    cout << "\nPress Enter to exit...";
    cin.ignore();
    cin.get();

    return 0;
}