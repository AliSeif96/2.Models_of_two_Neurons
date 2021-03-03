/************************************************************************************************/
/*** Topic: Two-Cell PING model with The E-cell is an RTM neuron, and the I-cell a WB neuron  ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 3/2/2021                                                                Ali-Seif   ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
//##############################################################
//####                                                      ####
//####               Calculate alpha and betas              ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________E-cell____________________________________________________
//________________________________________________________________________________________________
double alpha_n_e(double v_e) {
    return   0.032 * (v_e + 52.0) / (1.0 - exp(-(v_e + 52.0) / 5.0));
}
double beta_n_e(double v_e) {
    return   0.5 * exp(-(v_e + 57.0) / 40.0);
}
double alpha_m_e(double v_e) {
    return  0.32 * (v_e + 54.0) / (1.0 - exp(-(v_e + 54.0) / 4.0));
}
double beta_m_e(double v_e) {
    return  0.28 * (v_e + 27.0) / (exp((v_e + 27.0) / 5.0) - 1.0);
}
double alpha_h_e(double v_e) {
    return   0.128 * exp(-(v_e + 50) / 18);
}
double beta_h_e(double v_e) {
    return   4.0 / (1.0 + exp(-(v_e + 27.0) / 5.0));
}
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double alpha_n_i(double v_i) {
    return   -0.01 * (v_i + 34.0) / (exp(-0.1 * (v_i + 34.0)) - 1.0);
}
double beta_n_i(double v_i) {
    return   0.125 * exp(-(v_i + 44.0) / 80.0);
}
double alpha_m_i(double v_i) {
    return  0.1 * (v_i + 35.0) / (1.0 - exp(-(v_i + 35.0) / 10.0));
}
double beta_m_i(double v_i) {
    return  4.0 * exp(-(v_i + 60.0) / 18.0);
}
double alpha_h_i(double v_i) {
    return   0.07 * exp(-(v_i + 58.0) / 20.0);
}
double beta_h_i(double v_i) {
    return   1.0 / (exp(-0.1 * (v_i + 28.0)) + 1.0);
}
//##############################################################
//####                                                      ####
//####      Calculate infinite activation variables         ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________E-cell____________________________________________________
//________________________________________________________________________________________________
double n_inf_e(double v_e) {
    return   alpha_n_e(v_e) / (alpha_n_e(v_e) + beta_n_e(v_e));
}
double h_inf_e(double v_e) {
    return   alpha_h_e(v_e) / (alpha_h_e(v_e) + beta_h_e(v_e));
}
double m_inf_e(double v_e) {
    return   alpha_m_e(v_e) / (alpha_m_e(v_e) + beta_m_e(v_e));
}
double tau_n_e(double v_e) {
    return   (1.0 / (alpha_n_e(v_e) + beta_n_e(v_e)));
}
double tau_h_e(double v_e) {
    return   (1.0 / (alpha_h_e(v_e) + beta_h_e(v_e)));
}
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double n_inf_i(double v_i) {
    return   alpha_n_i(v_i) / (alpha_n_i(v_i) + beta_n_i(v_i));
}
double h_inf_i(double v_i) {
    return   alpha_h_i(v_i) / (alpha_h_i(v_i) + beta_h_i(v_i));
}
double m_inf_i(double v_i) {
    return   alpha_m_i(v_i) / (alpha_m_i(v_i) + beta_m_i(v_i));
}
double tau_n_i(double v_i) {
    double tau_n = 1.0 / (alpha_n_i(v_i) + beta_n_i(v_i));
    int phi = 5;
    return   (tau_n / phi);
}
double tau_h_i(double v_i) {
    double tau_h = 1.0 / (alpha_h_i(v_i) + beta_h_i(v_i));
    int phi = 5;
    return    (tau_h / phi);
}
//##############################################################
//####                                                      ####
//####             Calculation of currents                  ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________E-cell____________________________________________________
//________________________________________________________________________________________________
double INa_e(double v_e, double h_e) {
    int   g_Na_e = 100;
    int   E_Na_e = 50;
    return    g_Na_e * h_e * pow (m_inf_e(v_e), 3) * (v_e - E_Na_e);
}
double IK_e(double v_e, double n_e) {
    int   g_K = 80;
    int   E_K = -100;
    return    g_K * pow(n_e, 4) * (v_e - E_K);
}
double Il_e(double v_e) {
    float g_l = 0.1;
    int   E_l = -67;
    return        g_l * (v_e - E_l);
}
double Isyn_e(double v_e,double g_ie,double s_i) {
    int   E_rev_i = -75;
    return        g_ie * s_i * (E_rev_i - v_e);
}
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double INa_i(double v_i, double h_i) {
    int   g_Na_i = 35;
    int   E_Na_i = 55;
    return    g_Na_i * h_i * pow(m_inf_i(v_i), 3) * (v_i - E_Na_i);
}
double IK_i(double v_i, double n_i) {
    int   g_K = 9;
    int   E_K = -90;
    return    g_K * pow(n_i, 4) * (v_i - E_K);
}
double Il_i(double v_i) {
    float g_l = 0.1;
    int   E_l = -65;
    return        g_l * (v_i - E_l);
}
double Isyn_i(double v_i, double g_ei, double s_e) {
    int   E_rev_e = 0;
    return        g_ei * s_e * (E_rev_e - v_i);
}

//##############################################################
//####                                                      ####
//####            Differential Equations                    ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________E-cell____________________________________________________
//________________________________________________________________________________________________
double dvdt_e(double t, double v_e, double n_e, double h_e,double g_ie,double s_i) {
    float I_app_e = 1.4;
    float C_m = 1.0;
    return   (1 / C_m) * (I_app_e - (INa_e( v_e, h_e) + IK_e(v_e, n_e) + Il_e(v_e)) + Isyn_e(v_e,g_ie,s_i));
}
double dndt_e(double t, double n_e, double v_e) {
    return   ((n_inf_e(v_e)  - n_e)/ tau_n_e(v_e));
}
double dhdt_e(double t, double h_e, double v_e) {
    return   ((h_inf_e(v_e) - h_e) / tau_h_e(v_e));
}
double dqdt_e(double t, double q_e, double v_e,double tau_dq_e) {
    return   0.5 * (1 + tanh(0.1 * v_e)) * (1.0 - q_e) * 10.0 - q_e / tau_dq_e;
}
double dsdt_e(double t, double s_e,double q_e, double v_e, double tau_r_e, double tau_d_e) {
    return   q_e * (1.0 - s_e) / tau_r_e - s_e / tau_d_e;
}
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double dvdt_i(double t, double v_i, double n_i, double h_i, double g_ei, double s_e) {
    float I_app_i = 0;
    float C_m = 1.0;
    return   (1 / C_m) * (I_app_i - (INa_i(v_i, h_i) + IK_i(v_i, n_i) + Il_i(v_i)) + Isyn_i(v_i, g_ei, s_e));
}
double dndt_i(double t, double n_i, double v_i) {
    return   ((n_inf_i(v_i) - n_i) / tau_n_i(v_i));
}
double dhdt_i(double t, double h_i, double v_i) {
    return   ((h_inf_i(v_i) - h_i) / tau_h_i(v_i));
}
double dqdt_i(double t, double q_i, double v_i, double tau_dq_i) {
    return   0.5 * (1 + tanh(0.1 * v_i)) * (1.0 - q_i) * 10.0 - q_i / tau_dq_i;
}
double dsdt_i(double t, double s_i, double q_i, double v_i, double tau_r_i, double tau_d_i) {
    return   q_i * (1.0 - s_i) / tau_r_i - s_i / tau_d_i;
}

//##############################################################
//####                                                      ####
//####               tau_dq_function                        ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//_____________________________________tau_peak___________________________________________________
//________________________________________________________________________________________________
double tau_peak_function(double tau_d, double tau_r, double tau_d_q) {

    double dt = 0.01;
    double dt05 = dt / 2;

    double s = 0;
    double t = 0;
    double s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s * tau_d;

    double t_old;
    double s_inc_old;
    double s_tmp;
    double s_inc_tmp;
    while (s_inc > 0) {
        t_old = t;
        s_inc_old = s_inc;
        s_tmp = s + dt05 * s_inc;
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * (1 - s_tmp) / tau_r - s_tmp / tau_d;
        s = s + dt * s_inc_tmp;
        t = t + dt;
        s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s / tau_d;
    }
    return   (t_old * (-s_inc) + t * s_inc_old) / (s_inc_old - s_inc);
}
//________________________________________________________________________________________________
//______________________________________tau_dq____________________________________________________
//________________________________________________________________________________________________
double tau_dq_function(double tau_d, double tau_r, double tau_hat) {
    double tau_d_q_left = 1.0;
    while (tau_peak_function( tau_d,  tau_r, tau_d_q_left) > tau_hat) {
        tau_d_q_left = tau_d_q_left / 2;
    }

    double tau_d_q_right = tau_r;
    while (tau_peak_function(tau_d, tau_r, tau_d_q_right) < tau_hat) {
        tau_d_q_right = tau_d_q_right * 2;
    }
    double tau_d_q_mid;

    while (tau_d_q_right - tau_d_q_left > pow(10, -12)) {
        tau_d_q_mid = (tau_d_q_left + tau_d_q_right) / 2;
        if (tau_peak_function(tau_d, tau_r, tau_d_q_mid) <= tau_hat) {
            tau_d_q_left = tau_d_q_mid;

        }
        else {
            tau_d_q_right = tau_d_q_mid;

        }
    }
    return (tau_d_q_left + tau_d_q_right) / 2;
}


//_______________________________________________________________________________________\\
//_____________              The principle of the program                   _____________\\
//_____________                                      @                      _____________\\
//_____________           @@       @@       @            @@     @           _____________\\
//_____________           @ @     @ @      @ @       @   @ @    @           _____________\\
//_____________           @  @   @  @     @   @      @   @  @   @           _____________\\
//_____________           @   @@@   @    @@@@@@@     @   @   @  @           _____________\\
//_____________           @    @    @   @       @    @   @    @ @           _____________\\
//_____________           @         @  @         @   @   @     @@           _____________\\
//_______________________________________________________________________________________
int main() {

    double  t0 = 0, t_final = 201, dt = 0.01, dt2 = dt / 2;
    //________________________________________________________________________________________________
    //______________________________________E-cell____________________________________________________
    //________________________________________________________________________________________________
    double  I_e = 1.4,   
            g_ei = 0.25, 
            v_rev_e = 0, 
            tau_r_e = 0.5, 
            tau_peak_e = 0.5,
            tau_d_e = 3;
    double  tau_dq_e = tau_dq_function(tau_d_e, tau_r_e, tau_peak_e);//0.1723;
    //________________________________________________________________________________________________
    //______________________________________I-cell____________________________________________________
    //________________________________________________________________________________________________
    double  I_i = 0, 
            g_ie = 0.25, 
            v_rev_i = -75.0, 
            tau_r_i = 0.5, 
            tau_peak_i = 0.5,
            tau_d_i = 9;
    double  tau_dq_i = tau_dq_function(tau_d_i, tau_r_i, tau_peak_i); //0.1163;
    //initial conditions
    double  v_e = -75.0,
            h_e = 0.1,
            n_e = 0.1,
            q_e = 0,
            s_e = 0,
            v_i = -75.0,
            h_i = 0.1,
            n_i = 0.1,
            q_i = 0,
            s_i = 0;

    ofstream temp("temp.txt", ios::out | ios::trunc);

    double v_inc_e = 0.0, v_temp_e = 0.0, h_inc_e = 0.0, h_temp_e = 0.0, n_inc_e = 0.0, n_temp_e = 0.0, q_inc_e = 0.0, q_temp_e = 0.0, s_inc_e = 0.0, s_temp_e = 0.0;
    double v_inc_i = 0.0, v_temp_i = 0.0, h_inc_i = 0.0, h_temp_i = 0.0, n_inc_i = 0.0, n_temp_i = 0.0, q_inc_i = 0.0, q_temp_i = 0.0, s_inc_i = 0.0, s_temp_i = 0.0;

    for (t0 = 0; t0 <= t_final; t0 = t0 + dt) {

        v_inc_e = dvdt_e(t0, v_e, n_e, h_e,g_ie,s_i);
        n_inc_e = dndt_e(t0, n_e, v_e);
        h_inc_e = dhdt_e(t0, h_e, v_e);
        q_inc_e = dqdt_e(t0, q_e, v_e,tau_dq_e);
        s_inc_e = dsdt_e(t0,s_e,q_e,v_e,tau_r_e,tau_d_e);

        v_inc_i = dvdt_i(t0, v_i, n_i, h_i, g_ei, s_e);
        n_inc_i = dndt_i(t0, n_i, v_i);
        h_inc_i = dhdt_i(t0, h_i, v_i);
        q_inc_i = dqdt_i(t0, q_i, v_i, tau_dq_i);
        s_inc_i = dsdt_i(t0, s_i, q_i, v_i, tau_r_i, tau_d_i);


        v_temp_e = v_e + dt2 * v_inc_e;
        h_temp_e = h_e + dt2 * h_inc_e;
        n_temp_e = n_e + dt2 * n_inc_e;
        q_temp_e = q_e + dt2 * q_inc_e;
        s_temp_e = s_e + dt2 * s_inc_e;


        v_temp_i = v_i + dt2 * v_inc_i;
        h_temp_i = h_i + dt2 * h_inc_i;
        n_temp_i = n_i + dt2 * n_inc_i;
        q_temp_i = q_i + dt2 * q_inc_i;
        s_temp_i = s_i + dt2 * s_inc_i;


        v_inc_e = dvdt_e(t0, v_temp_e, n_temp_e, h_temp_e, g_ie, s_temp_i);
        n_inc_e = dndt_e(t0, n_temp_e, v_temp_e);
        h_inc_e = dhdt_e(t0, h_temp_e, v_temp_e);
        q_inc_e = dqdt_e(t0, q_temp_e, v_temp_e, tau_dq_e);
        s_inc_e = dsdt_e(t0, s_temp_e, q_temp_e, v_temp_e, tau_r_e, tau_d_e);

        v_inc_i = dvdt_i(t0, v_temp_i, n_temp_i, h_temp_i, g_ei, s_temp_e);
        n_inc_i = dndt_i(t0, n_temp_i, v_temp_i);
        h_inc_i = dhdt_i(t0, h_temp_i, v_temp_i);
        q_inc_i = dqdt_i(t0, q_temp_i, v_temp_i, tau_dq_i);
        s_inc_i = dsdt_i(t0, s_temp_i, q_temp_i, v_temp_i, tau_r_i, tau_d_i);



        v_e = v_e + dt * v_inc_e;
        h_e = h_e + dt * h_inc_e;
        n_e = n_e + dt * n_inc_e;
        q_e = q_e + dt * q_inc_e;
        s_e = s_e + dt * s_inc_e;

        v_i = v_i + dt * v_inc_i;
        h_i = h_i + dt * h_inc_i;
        n_i = n_i + dt * n_inc_i;
        q_i = q_i + dt * q_inc_i;
        s_i = s_i + dt * s_inc_i;

        temp << t0 << '\t' << v_e << '\t' << v_i << endl;
    }
    temp.close();
    cout << "\nFinish" << endl;
    return 0;
}