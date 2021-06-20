






//                                C++ CODE FOR THE NUMERICAL CALCULATION OF STELLAR INTERIORS
//
// The code below is derived from the equations in Papers 2 and 3 of the series :
//
// Paper 1 : TechnicalNote_StellarStructure         (this is TBD, not uploaded yet)
// Paper 2 : TechnicalNote_NumericalIntegration     (this is uploaded and substantially complete)
// Paper 3 : TechnicalNote_CplusplusCode            (this is uploaded but lacks sections on transformed variables and Henyey calculations)
//
//
// The code below does essentially four things :
// Calculates a few starting values at the stellar centre using a modified constant density approximation
// Inputs these starting values into an Adams-Bashforth integration algorithm and calculates outwards towards the surface of the star
// Calculates a few starting values at the stellar surface using an approximation where mass and luminosity are assumed constant and at their full values
// Inputs these starting values into the same Adams-Bashforth integration algorithm and calculates inwards towards the centre of the star
// There is some undeveloped code for calculations in transformed variables. This needs completion.
// The code is C++ and was developed under VisualStudio v16.8.4
//
// Figure 1 in Paper 3 gives the basic differential equations
// Figure 2 in Paper 3 gives the constant density approximation equations for the stellar centre
// Section 2.1.2 in Paper 3 gives the constant mass and luminosity approximation equations for the stellar surface
// Section 2.1.3 in Paper 3 describes the Adams-Bashforth integration algorithm 
//
// The code below does not :
// Attempt to join the outward and inward solutions to get a full solution
// Do calculations for convective energy transfer
// Do evolutionary calculations where the star develops in time and its chemical composition changes
// Apply the Henyey method of stellar structure calculation
// Do calculations with the stellar mass (instead of distance) as the independent variable.
//
// This code is not guaranteed to be bug-free nor to be in the most elegant form.
// Please use and develop the uploaded Papers 1-3 and this code as you wish but credit me please : 
//
// Peter Gardner : peter.gardner13@btopenworld.com, UK 07495 416625, also see my LinkedIn profile.
//
// 

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#define calculation_length 190

const int numberofsteps = 10;
const int ab_calculation_length = calculation_length + numberofsteps;
double ab_p_values_out[ab_calculation_length] = { 0.0 };
double ab_p_derivatives_out[ab_calculation_length] = { 0.0 };
double ab_q_values_out[ab_calculation_length] = { 0.0 };
double ab_q_derivatives_out[ab_calculation_length] = { 0.0 };
double ab_f_values_out[ab_calculation_length] = { 0.0 };
double ab_f_derivatives_out[ab_calculation_length] = { 0.0 };
double ab_t_values_out[ab_calculation_length] = { 0.0 };
double ab_t_derivatives_out[ab_calculation_length] = { 0.0 };

double ab_p_values_in[ab_calculation_length] = { 0.0 };
double ab_p_derivatives_in[ab_calculation_length] = { 0.0 };
double ab_q_values_in[ab_calculation_length] = { 0.0 };
double ab_q_derivatives_in[ab_calculation_length] = { 0.0 };
double ab_f_values_in[ab_calculation_length] = { 0.0 };
double ab_f_derivatives_in[ab_calculation_length] = { 0.0 };
double ab_t_values_in[ab_calculation_length] = { 0.0 };
double ab_t_derivatives_in[ab_calculation_length] = { 0.0 };


int main(int argc, char* argv[]) {


    ofstream OutResults_File_pctc("Out_Results_File_pctc.txt");
    ofstream OutResults_File_mrl("Out_Results_File_mrl.txt");


    //FUNDAMENTAL CONSTANTS
    //CGS UNITS (Centimetre, Gram, Second)

    const double gravitational_constant = 6.67 * pow(10, -8);
    const double speed_of_light = 3.0 * pow(10, 10);
    const double radiation_density_constant = 7.56 * pow(10, -15);      //"a". Schwarszchild page 41. The book apparently mistakenly calls it the Stefan Boltzmann constant.
    const double hydrogen_fraction = 0.98;
    const double helium_fraction = 0.02;
    const double heavy_fraction = 0.00;
    const double boltzmann_constant = 1.38 * pow(10, -16);
    const double proton_mass = 1.67 * pow(10, -24);
    const double pi = 3.14159;

    //DERIVED CONSTANTS

    double mu_reciprocal = 2.0 * hydrogen_fraction + 0.75 * helium_fraction + 0.5 * heavy_fraction;       //This is fixed by assumption.
    double mu = 1.0 / mu_reciprocal;
    const double t_over_tentosix = 12.0;                          //Schwarzschild page 83, Table 10.1
    const double epsilon1 = pow(10, -5.5);                        //Schwarzschild page 83, Table 10.1
    const double nu = 4.5;                                        //Schwarzschild page 83, Table 10.1

    double epp_tentosix_over_t = 0.0;
    double epp_tentosix_over_t_to_twothirds = 0.0;
    double epp_pressure_term = 0.0;
    double epp_exponential_power = 0.0;
    double epp_exponential_term = 0.0;
    double energy_generation = 0.0;

    const double gbart = 1.0;                                                                             // Schwarzschild equation 9.16
    const double gffbar = 1.0;                                                                            // Schwarzschild equation 9.18
    bool gas_equation_density = true;
    double central_density = 0.0;
    double chi_zero = 0.0;
    double chi_bf = 0.0;
    double chi_ff = 0.0;
    double central_opacity = 0.0;


    const int numberofsteps = 10;
    bool adams_bashforth = true;
    bool untransformed = true;
    bool transformed = false;
    bool linear_density_decrease = true;
    double step = 0.0;

    if (untransformed) {
        step = 3.5 * pow(10, 8);
    }
    else {                                                 //step size for untransformed variables
        step = 0.1;                                        //step size for transformed variables
    }

    double new_density = 0.0;
    double density = 0.0;

    double lhs = 0.0;
    double rhs = 0.0;
    double r = 0.0;
    double tmp, tmpr;
    const bool pctc = true;
    const bool mrl = true;

    chi_bf = (4.34 * pow(10, 25) * gbart * heavy_fraction * (1.0 + hydrogen_fraction));                              // Schwarzschild equation 9.16
    chi_ff = (3.68 * pow(10, 22) * gffbar * (hydrogen_fraction + helium_fraction) * (1.0 + hydrogen_fraction));      // Schwarzschild equation 9.18
    chi_zero = chi_bf + chi_ff;

    if (pctc) {

        // PRESSURE LOOP

        const int numberofsteps_pressure = 5;
        const int numberofsteps_temperature = 5;
        //Schwarszchild equation 8.2
        double Central_Pressure = 0.0;
        double Central_Temperature = 0.0;
        double pressure_exponent = 0.0;
        double pressure_multiplier = 2.0;
        double temperature_exponent = 0.0;
        double temperature_multiplier = 0.0;

        for (int m = 0; m < numberofsteps_pressure; m++) {

            pressure_exponent = 17.0;
            pressure_multiplier = double(m) + 1.0;
            Central_Pressure = pressure_multiplier * pow(10, pressure_exponent);

            for (int k = 0; k < numberofsteps_temperature; k++) {

                temperature_exponent = 7.0;
                temperature_multiplier = double(k) + 2.0;
                Central_Temperature = temperature_multiplier * pow(10, temperature_exponent);

                OutResults_File_pctc << Central_Pressure;
                OutResults_File_pctc << "          ";
                OutResults_File_pctc << Central_Temperature;
                OutResults_File_pctc << "\n";

                std::cout << "Pressure_Exponent =  " << pressure_exponent;
                std::cout << "\n";
                std::cout << "Central Pressure =  " << Central_Pressure;
                std::cout << "\n";

                if (gas_equation_density == true) {
                    std::cout << "\n";
                    central_density = (mu * Central_Pressure * proton_mass) / (boltzmann_constant * Central_Temperature);
                    std::cout << "\n";
                    std::cout << "\n";
                }


                if (gas_equation_density == false) {
                    std::cout << "\n";
                    central_density = 72.0;
                    std::cout << "\n";
                    std::cout << "\n";
                }

                central_opacity = chi_zero * central_density / pow(Central_Temperature, 3.5);

                double central_energy_generation = 0.0;
                const double interpolated_central_energy_generation = epsilon1 * central_density * hydrogen_fraction * hydrogen_fraction * pow(t_over_tentosix, nu);


                std::cout << "Composition =  " << mu;
                std::cout << "\n";
                std::cout << "Central Pressure =  " << Central_Pressure;
                std::cout << "\n";
                std::cout << "Central Temperature =  " << Central_Temperature;
                std::cout << "\n";
                std::cout << "Central Density =  " << central_density;
                std::cout << "\n";
                std::cout << "Interpolated_Central Energy Generation =  " << interpolated_central_energy_generation;
                std::cout << "\n";
                std::cout << "Chi_bf =  " << chi_bf;
                std::cout << "\n";
                std::cout << "Chi_ff =  " << chi_ff;
                std::cout << "\n";
                std::cout << "Central Opacity =  " << central_opacity;
                std::cout << "\n";
                std::cout << "\n";

                double pressure, mass, luminosity, temperature;

                // VARIABLES OF THE ALGORITHM (OUTWARDS SOLUTION)

                double p_centre[numberofsteps];
                double q_centre[numberofsteps];
                double f_centre[numberofsteps];
                double t_centre[numberofsteps];

                double p_centre_derivatives[numberofsteps];
                double q_centre_derivatives[numberofsteps];
                double f_centre_derivatives[numberofsteps];
                double t_centre_derivatives[numberofsteps];

                // VARIABLES OF THE ALGORITHM (TRANSFORMED VARIABLES)

                double xstar, pstar, qstar, fstar, tstar;

                // Constant C : Schwarszchild equation 13.9
                //double c1 = 3.0 / (4.0 * radiation_density_constant * speed_of_light);
                //double c2 = boltzmann_constant / (proton_mass * gravitational_constant);
                //double c3 = pow(c2, 7.5);
                //double c4 = (1.0) / (pow(4.0 * pi, 3));
                //double c5 = central_opacity / (pow(mu, 7.5));
                //double c6 = Model_Luminosity * pow(Model_Radius, 0.5) / pow(Model_Mass, 5.5);
                //double C = c1 * c3 * c4 * c5 * c6;

                // Constant D : Schwarszchild equation 13.10
                //double d1 = proton_mass * gravitational_constant / boltzmann_constant;
                //double d2 = pow(d1, nu) / (4.0 * pi);
                //double d3 = central_energy_generation * pow(mu, nu);
                //double d4 = pow(Model_Mass, nu + 2) / (Model_Luminosity * pow(Model_Radius, nu + 3));
                //double D = d2 * d3 * d4;

                // The starred variables : Schwarszchild equations 13.13 
                //double tc = (boltzmann_constant * Model_Radius * Central_Temperature) / (mu * proton_mass * gravitational_constant * Model_Mass);
                //double p_zero_cubed = pow(tc, 9.5 - nu) / (C * D);
                //double p_zero = pow(p_zero_cubed, 1.0/3.0);
                //double p_zero = cbrt(p_zero_cubed);
                //double pcentre_numerator = 4.0 * pi * pow(Model_Radius, 4);
                //double pcentre_denominator = gravitational_constant * pow(Model_Mass, 2);
                //double pcentre = (pcentre_numerator / pcentre_denominator) * Central_Pressure;
                //The value of pstar at the centre
                //double pstarc = pcentre / p_zero;

                //double x_zero_1 = pow(Central_Temperature, 9.5 - nu + 2) / (C * D);
                //double x_zero = pow(x_zero_1, 0.5) / pow(p_zero, 2);

                //Atomic abundance profile (H, He, heavies)
                //const int NumberOfLayers = 4;
                //const double Layer_Width = Model_Radius / double(NumberOfLayers);
                //double AH[NumberOfLayers] = { 0.08, 0.53, 0.20, 0.06 };
                //double BH[NumberOfLayers] = { 0.34, 0.40, 0.80, 0.95 };           // apply factor of pow(10, -10)
                //double AHe[NumberOfLayers] = { 0.08, 0.50, 0.20, 0.06 };
                //double BHe[NumberOfLayers] = { -0.64, -0.58, -0.20, -0.05 };     // apply factor of pow(10, -10)
                //double H_abund = 0.0, He_abund = 0.0, Hv_abund = 0.0;

                //THE SUN

                //mass = 2 * pow(10, 33) grams
                //radius = 7 * pow(10, 10) cms
                //Central_Pressure = 6.0 * pow(10, 15)
                //central_density = 162 grams/cm3
                //Central_Temperature = 1.57 * pow(10, 7) Kelvin
                //luminosity = 3.8 * pow(10, 26) watts

                if (gas_equation_density == true) {
                    std::cout << "\n";
                    std::cout << "IDEAL GAS EQUATION CENTRAL DENSITY";
                    std::cout << "\n";
                    std::cout << "\n";
                }

                if (gas_equation_density == false) {
                    std::cout << "\n";
                    std::cout << "INTERMEDIATE CENTRAL DENSITY";
                    std::cout << "\n";
                    std::cout << "\n";
                }

                if (linear_density_decrease == true) {
                    std::cout << "LINEAR DENSITY DECREASE";
                    std::cout << "TWO TIMES STEP SIZE";
                }
                if (linear_density_decrease == false) {
                    std::cout << "CONSTANT DENSITY APPROXIMATION";
                }


                std::cout << "\n";
                std::cout << "\n";
                std::cout << "\n";

                if (untransformed) {

                    std::cout << "Untransformed variables";

                    std::cout << "\n";
                    std::cout << "\n";
                    std::cout << "\n";


                    // Schwarzschild page 81
                    epp_tentosix_over_t = 1000000.0 / Central_Temperature;
                    epp_tentosix_over_t_to_twothirds = pow(epp_tentosix_over_t, 2.0 / 3.0);
                    epp_pressure_term = 2.5 * pow(10, 6) * central_density * hydrogen_fraction * hydrogen_fraction;
                    epp_exponential_power = -33.8 * pow(epp_tentosix_over_t, 1.0 / 3.0);
                    epp_exponential_term = exp(epp_exponential_power);
                    central_energy_generation = epp_pressure_term * epp_tentosix_over_t_to_twothirds * epp_exponential_term;
                    tmp = central_energy_generation;
                    std::cout << setw(20) << tmp;

                    std::cout << "\n";
                    std::cout << "\n";
                    std::cout << "\n";

                    for (int n = 0; n < numberofsteps; n++) {
                        r = double(n) * step;

                        if (gas_equation_density == true) {
                            if (n == 0) {
                                // density at the centre according to the ideal gas law
                                density = (mu * Central_Pressure * proton_mass) / (boltzmann_constant * Central_Temperature);
                            }

                            if (n > 0) {
                                // density at the new distance according to the ideal gas law
                                density = (mu * pressure * proton_mass) / (boltzmann_constant * temperature);
                            }
                        }

                        if (gas_equation_density == false) {
                            density = 72.0;
                        }

                        pressure = Central_Pressure - (2.0 * pi * gravitational_constant * pow(density, 2) * pow(r, 2)) / 3.0;
                        mass = (4.0 * pi * density * pow(r, 3)) / 3.0;
                        luminosity = (4.0 * pi * density * central_energy_generation * pow(r, 3)) / 3.0;
                        temperature = Central_Temperature - (central_opacity * central_energy_generation * pow(density, 2) * pow(r, 2)) / (8.0 * radiation_density_constant * speed_of_light * pow(Central_Temperature, 3));

                        p_centre[n] = pressure;
                        q_centre[n] = mass;
                        f_centre[n] = luminosity;
                        t_centre[n] = temperature;

                        if (n == 0) {
                            p_centre_derivatives[0] = 0.0;
                            q_centre_derivatives[0] = 0.0;
                            f_centre_derivatives[0] = 0.0;
                            t_centre_derivatives[0] = 0.0;
                        }

                        if (n > 0) {
                            p_centre_derivatives[n] = -(density * gravitational_constant * q_centre[n]) / pow(r, 2);
                            q_centre_derivatives[n] = density * 4.0 * pi * pow(r, 2);
                            f_centre_derivatives[n] = density * 4.0 * pi * pow(r, 2) * central_energy_generation;
                            t_centre_derivatives[n] = -(3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (central_opacity * density / pow(t_centre[n], 3)) * (f_centre[n] / (4.0 * pi * pow(r, 2)));
                        }

                        std::cout << setw(10) << r << "    ";
                        std::cout << setw(5) << density << "    ";
                        std::cout << setw(5) << p_centre[n] << "    ";
                        std::cout << setw(5) << q_centre[n] << "    ";
                        std::cout << setw(5) << f_centre[n] << "    ";
                        std::cout << setw(5) << t_centre[n] << "    ";
                        std::cout << setw(5) << p_centre_derivatives[n] << "    ";
                        std::cout << setw(5) << q_centre_derivatives[n] << "    ";
                        std::cout << setw(5) << f_centre_derivatives[n] << "    ";
                        std::cout << setw(10) << t_centre_derivatives[n] << "    ";
                        std::cout << "\n";

                    }
                }

                if (transformed) {

                    std::cout << "Transformed variables";

                    std::cout << "\n";
                    std::cout << "\n";
                    std::cout << "\n";

                    for (int n = 0; n < numberofsteps; n++) {
                        xstar = double(n) * step;

                        //pstar = pstarc - pow(pstarc, 2) * pow(xstar, 2) / 6.0;
                        //qstar = pstarc * pow(xstar, 3) / 3.0;
                        //fstar = pow(pstarc, 2) * pow(xstar, 3) / 3.0;
                        //tstar = 1.0 - pow(pstarc, 4) * pow(xstar, 2) / 6.0;

                        //p_centre[n] = pstar;
                        //q_centre[n] = qstar;
                        //f_centre[n] = fstar;
                        //t_centre[n] = tstar;

                        //std::cout << setw(20) << xstar << "    ";
                        //std::cout << setw(20) << pstar << "    ";
                        //std::cout << setw(20) << qstar << "    ";
                        //std::cout << setw(20) << fstar << "    ";
                        //std::cout << setw(20) << tstar << "    ";
                        //std::cout << "\n";
                    }
                }


                std::cout << "\n";
                std::cout << "\n";
                std::cout << "\n";

                if (adams_bashforth) {
                    std::cout << "\n";
                    std::cout << "\n";
                    std::cout << "\n";
                    std::cout << "\n";
                    std::cout << "INTEGRATING OUTWARD USING ADAMS-BASHFORTH METHOD";
                    std::cout << "\n";

                    for (int k = 0; k < numberofsteps; k++) {
                        ab_p_values_out[k] = p_centre[k];
                        ab_p_derivatives_out[k] = p_centre_derivatives[k];
                    }

                    for (int k = 0; k < numberofsteps; k++) {
                        ab_q_values_out[k] = q_centre[k];
                        ab_q_derivatives_out[k] = q_centre_derivatives[k];
                    }

                    for (int k = 0; k < numberofsteps; k++) {
                        ab_f_values_out[k] = f_centre[k];
                        ab_f_derivatives_out[k] = f_centre_derivatives[k];
                    }

                    for (int k = 0; k < numberofsteps; k++) {
                        ab_t_values_out[k] = t_centre[k];
                        ab_t_derivatives_out[k] = t_centre_derivatives[k];
                    }

                    OutResults_File_pctc << setw(20) << "CENTRAL PRESSURE = " << Central_Pressure;
                    OutResults_File_pctc << setw(20);
                    OutResults_File_pctc << setw(20) << "CENTRAL TEMPERATURE = " << Central_Temperature;
                    OutResults_File_pctc << "\n";
                    OutResults_File_pctc << "\n";


                    for (int k = numberofsteps; k < ab_calculation_length / 2; k++) {

                        tmpr = double(k) * step;

                        //Adams-Bashforth third order predictor method

                        ab_p_values_out[k] = ab_p_values_out[k - 1] + step * (23.0 * ab_p_derivatives_out[k - 1] - 16.0 * ab_p_derivatives_out[k - 2] + 5.0 * ab_p_derivatives_out[k - 3]) / 12.0;
                        if (ab_p_values_out[k] < 0.0) {
                            std::cout << "\n";
                            std::cout << setw(20) << "PRESSURE BELOW ZERO";
                            OutResults_File_pctc << "PRESSURE BELOW ZERO";
                            std::cout << "\n";
                            break;
                        }

                        std::cout << setw(20) << tmpr;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmpr;
                        OutResults_File_pctc << "P";

                        tmp = ab_p_values_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;
                        OutResults_File_pctc << "P";

                        tmp = ab_p_derivatives_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_p_derivatives_out[k - 2];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_p_derivatives_out[k - 3];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_p_values_out[k];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        std::cout << "\n";
                        OutResults_File_pctc << "\n";

                        ab_q_values_out[k] = ab_q_values_out[k - 1] + step * (23.0 * ab_q_derivatives_out[k - 1] - 16.0 * ab_q_derivatives_out[k - 2] + 5.0 * ab_q_derivatives_out[k - 3]) / 12.0;

                        if (ab_q_values_out[k] < 0.0) {
                            std::cout << "\n";
                            std::cout << setw(20) << "MASS BELOW ZERO";
                            OutResults_File_pctc << "MASS BELOW ZERO";
                            std::cout << "\n";
                            break;
                        }

                        std::cout << setw(20) << tmpr;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmpr;
                        OutResults_File_pctc << "Q";

                        tmp = ab_q_values_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;
                        OutResults_File_pctc << "Q";

                        tmp = ab_q_derivatives_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_q_derivatives_out[k - 2];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_q_derivatives_out[k - 3];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_q_values_out[k];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        std::cout << "\n";
                        OutResults_File_pctc << "\n";


                        ab_f_values_out[k] = ab_f_values_out[k - 1] + step * (23.0 * ab_f_derivatives_out[k - 1] - 16.0 * ab_f_derivatives_out[k - 2] + 5.0 * ab_f_derivatives_out[k - 3]) / 12.0;

                        if (ab_f_values_out[k] < 0.0) {
                            std::cout << "\n";
                            std::cout << setw(20) << "LUMINOSITY BELOW ZERO";
                            OutResults_File_mrl << "LUMINOSITY BELOW ZERO";
                            std::cout << "\n";
                            break;
                        }

                        std::cout << setw(20) << tmpr;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmpr;
                        OutResults_File_pctc << "F";

                        tmp = ab_f_values_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;
                        OutResults_File_pctc << "F";

                        tmp = ab_f_derivatives_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_f_derivatives_out[k - 2];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_f_derivatives_out[k - 3];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_f_values_out[k];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        std::cout << "\n";
                        OutResults_File_pctc << "\n";


                        ab_t_values_out[k] = ab_t_values_out[k - 1] + step * (23.0 * ab_t_derivatives_out[k - 1] - 16.0 * ab_t_derivatives_out[k - 2] + 5.0 * ab_t_derivatives_out[k - 3]) / 12.0;

                        if (ab_t_values_out[k] < 0.0) {
                            std::cout << "\n";
                            std::cout << setw(20) << "TEMPERATURE BELOW ZERO";
                            OutResults_File_pctc << "TEMPERATURE BELOW ZERO";
                            std::cout << "\n";
                            break;
                        }

                        std::cout << setw(20) << tmpr;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmpr;
                        OutResults_File_pctc << "T";

                        tmp = ab_t_values_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;
                        OutResults_File_pctc << "T";

                        tmp = ab_t_derivatives_out[k - 1];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_t_derivatives_out[k - 2];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_t_derivatives_out[k - 3];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        tmp = ab_t_values_out[k];
                        std::cout << setw(20) << tmp;
                        OutResults_File_pctc << setw(20);
                        OutResults_File_pctc << tmp;

                        std::cout << "\n";
                        OutResults_File_pctc << "\n";

                        density = (mu * ab_p_values_out[k] * proton_mass) / (boltzmann_constant * ab_t_values_out[k]);
                        tmp = density;
                        std::cout << setw(20) << tmp;

                        // Schwarzschild page 81
                        epp_tentosix_over_t = 1000000.0 / ab_t_values_out[k];
                        epp_tentosix_over_t_to_twothirds = pow(epp_tentosix_over_t, 2.0 / 3.0);
                        epp_pressure_term = 2.5 * pow(10, 6) * density * hydrogen_fraction * hydrogen_fraction;
                        epp_exponential_power = -33.8 * pow(epp_tentosix_over_t, 1.0 / 3.0);
                        epp_exponential_term = exp(epp_exponential_power);
                        energy_generation = epp_pressure_term * epp_tentosix_over_t_to_twothirds * epp_exponential_term;
                        tmp = energy_generation;
                        std::cout << setw(20) << tmp;


                        //if (k == (ab_calculation_length - 1)) {

                        //    OutResults_File << density;
                        //    OutResults_File << "          ";
                        //    OutResults_File << energy_generation;
                        //    OutResults_File << "\n";
                        //    OutResults_File << ab_p_values_out[k - 1];
                        //    OutResults_File << "\n";
                        //    OutResults_File << ab_q_values_out[k - 1];
                        //    OutResults_File << "\n";
                        //    OutResults_File << ab_f_values_out[k - 1];
                        //    OutResults_File << "\n";
                        //    OutResults_File << ab_t_values_out[k - 1];
                        //    OutResults_File << "\n";
                        //    OutResults_File << "\n";
                        //}

                        ab_p_derivatives_out[k] = -(density * gravitational_constant * ab_q_values_out[k]) / pow(tmpr, 2);
                        ab_q_derivatives_out[k] = density * 4.0 * pi * pow(tmpr, 2);
                        ab_f_derivatives_out[k] = density * 4.0 * pi * pow(tmpr, 2) * energy_generation;
                        ab_t_derivatives_out[k] = -(3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (central_opacity * density / pow(ab_t_values_out[k], 3)) * (ab_f_values_out[k] / (4.0 * pi * pow(tmpr, 2)));
                        std::cout << "\n";
                    }
                }
            }
        }
    }

    // INWARD SOLUTION

    if (mrl) {

        double from_surface = 0.0;
        double rfactor, temperature_constants_factor;
        double pressure_fundamental_factors, pressure_stellar_factors, pressure_all_factors;
        double pressure, mass, luminosity, temperature;

        double p_surface[numberofsteps];
        double q_surface[numberofsteps];
        double f_surface[numberofsteps];
        double t_surface[numberofsteps];

        double p_surface_derivatives[numberofsteps] = { 0.0 };
        double q_surface_derivatives[numberofsteps] = { 0.0 };
        double f_surface_derivatives[numberofsteps] = { 0.0 };
        double t_surface_derivatives[numberofsteps] = { 0.0 };

        //Surface Values

        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "SURFACE VALUES GOING IN";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";

        double Model_Mass = 2.0 * pow(10, 33);                  //Wikipedia (grams). This is fixed by assumption.
        double Model_Radius = 7.0 * pow(10, 10);
        double Model_Luminosity = 6.0 * pow(10, 33);            //Schwarzschild page 42. This is a parameter of the model

        if (untransformed) {

            std::cout << "UNTRANSFORMED VARIABLES";

            std::cout << "\n";
            std::cout << "\n";
            std::cout << "\n";

            temperature_constants_factor = (mu * proton_mass * gravitational_constant * Model_Mass) / (4.25 * boltzmann_constant);        //radiative
            //constants_factor = (mu * proton_mass * gravitational_constant * mass) / (2.5 * boltzmann_constant);       //convective
            pressure_fundamental_factors = (2.0 / 8.5) * (4.0 * radiation_density_constant * speed_of_light / 3.0) * (boltzmann_constant / proton_mass);
            pressure_stellar_factors = (4.0 * pi * gravitational_constant) / (mu * chi_zero);

            std::cout << setw(20) << mu << "    ";
            std::cout << setw(20) << chi_zero << "    ";
            std::cout << setw(20) << temperature_constants_factor << "    ";
            std::cout << setw(20) << pressure_fundamental_factors << "    ";
            std::cout << setw(20) << pressure_stellar_factors << "    ";
            std::cout << "\n";
            std::cout << "\n";

            double surface_density = 0.0;
            double surface_pressure = 0.0;
            double surface_temperature = 0.0;
            double surface_opacity = 0.0;
            //const double surface_t_over_tentosix = 5.0;                                                                      //Schwarzschild page 83, Table 10.1
            //const double surface_epsilon1 = 0.0000001;                                                                       //Schwarzschild page 83, Table 10.1
            //const double surface_nu = 6.0;                                                                                   //Schwarzschild page 83, Table 10.1
            double surface_energy_generation = 0.0;


            const int numberofsteps_radius = 3;
            double radius_multiplier = 0.0;
            const int numberofsteps_luminosity = 3;
            double luminosity_multiplier = 0.0;

            for (int j = 0; j < numberofsteps_radius; j++) {

                radius_multiplier = 6.0 + double(j);
                Model_Radius = radius_multiplier * pow(10, 10);

                for (int m = 0; m < numberofsteps_luminosity; m++) {

                    luminosity_multiplier = 5.0 + double(m);
                    Model_Luminosity = luminosity_multiplier * pow(10, 33);

                    std::cout << setw(20) << Model_Luminosity << "    ";
                    std::cout << "\n";
                    std::cout << "\n";

                    for (int n = 0; n < numberofsteps; n++) {
                        from_surface = double(n) * step;
                        r = Model_Radius - from_surface;
                        if (r < 0) {
                            break;
                        }
                        rfactor = (1.0 / r) - (1.0 / Model_Radius);

                        pressure_all_factors = pow((pressure_fundamental_factors * pressure_stellar_factors) * (Model_Mass / Model_Luminosity), 0.5);

                        surface_temperature = temperature_constants_factor * rfactor;
                        surface_pressure = pressure_all_factors * pow(surface_temperature, 4.25);                 //radiative

                        p_surface[n] = surface_pressure;
                        q_surface[n] = Model_Mass;
                        f_surface[n] = Model_Luminosity;
                        t_surface[n] = surface_temperature;

                        if (gas_equation_density == true) {
                            surface_density = (mu * surface_pressure * proton_mass) / (boltzmann_constant * surface_temperature);
                        }

                        if (gas_equation_density == false) {
                            // temporary estimate
                            surface_density = central_density / 60.0;
                        }

                        surface_opacity = chi_zero * surface_density / pow(surface_temperature, 3.5);

                        // Schwarzschild page 81
                        epp_tentosix_over_t = 1000000.0 / surface_temperature;
                        epp_tentosix_over_t_to_twothirds = pow(epp_tentosix_over_t, 2.0 / 3.0);
                        epp_pressure_term = 2.5 * pow(10, 6) * surface_density * hydrogen_fraction * hydrogen_fraction;
                        epp_exponential_power = -33.8 * pow(epp_tentosix_over_t, 1.0 / 3.0);
                        epp_exponential_term = exp(epp_exponential_power);
                        surface_energy_generation = epp_pressure_term * epp_tentosix_over_t_to_twothirds * epp_exponential_term;
                        tmp = surface_energy_generation;
                        //std::cout << setw(20) << tmp;
                        //std::cout << "\n";


                        if (n == 0) {
                            p_surface_derivatives[0] = 0.0;
                            q_surface_derivatives[0] = 0.0;
                            f_surface_derivatives[0] = 0.0;
                            t_surface_derivatives[0] = 0.0;
                        }

                        if (n > 0) {
                            // Note opposite signs from outward calculation.
                            p_surface_derivatives[n] = (surface_density * gravitational_constant * q_surface[n]) / pow(r, 2);
                            q_surface_derivatives[n] = -surface_density * 4.0 * pi * pow(r, 2);
                            f_surface_derivatives[n] = -surface_density * 4.0 * pi * pow(r, 2) * surface_energy_generation;
                            t_surface_derivatives[n] = (3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (surface_opacity * surface_density / pow(t_surface[n], 3)) * (f_surface[n] / (4.0 * pi * pow(r, 2)));
                        }


                        std::cout << setw(10) << r << " ";
                        std::cout << setw(5) << surface_density << "     ";
                        std::cout << setw(10) << surface_pressure << " ";
                        std::cout << setw(10) << surface_temperature << " ";
                        std::cout << setw(10) << p_surface_derivatives[n] << "  ";
                        std::cout << setw(10) << q_surface_derivatives[n] << "  ";
                        std::cout << setw(10) << f_surface_derivatives[n] << "  ";
                        std::cout << setw(10) << t_surface_derivatives[n] << "  ";
                        std::cout << "\n";
                        std::cout << "\n";
                    }

                    if (adams_bashforth) {
                        std::cout << "\n";
                        std::cout << "\n";
                        std::cout << "\n";
                        std::cout << "\n";
                        std::cout << "INTEGRATING INWARD USING ADAMS-BASHFORTH METHOD";
                        std::cout << "\n";


                        for (int k = 0; k < numberofsteps; k++) {
                            ab_p_values_in[k] = p_surface[k];
                            ab_p_derivatives_in[k] = p_surface_derivatives[k];
                        }

                        for (int k = 0; k < numberofsteps; k++) {
                            ab_q_values_in[k] = q_surface[k];
                            ab_q_derivatives_in[k] = q_surface_derivatives[k];
                        }

                        for (int k = 0; k < numberofsteps; k++) {
                            ab_f_values_in[k] = f_surface[k];
                            ab_f_derivatives_in[k] = f_surface_derivatives[k];
                        }

                        for (int k = 0; k < numberofsteps; k++) {
                            ab_t_values_in[k] = t_surface[k];
                            ab_t_derivatives_in[k] = t_surface_derivatives[k];
                        }


                        for (int k = numberofsteps; k < ab_calculation_length - 10; k++) {


                            //Adams-Bashforth third order predictor method

                            tmpr = Model_Radius - double(k) * step;
                            if (tmpr < 0.0) {
                                break;
                            }
                            //std::cout << setw(20) << tmpr;

                            ab_p_values_in[k] = ab_p_values_in[k - 1] + step * (23.0 * ab_p_derivatives_in[k - 1] - 16.0 * ab_p_derivatives_in[k - 2] + 5.0 * ab_p_derivatives_in[k - 3]) / 12.0;

                            if (ab_p_values_in[k] < 0.0) {
                                std::cout << "\n";
                                std::cout << setw(20) << "PRESSURE BELOW ZERO";
                                std::cout << "\n";
                                break;
                            }

                            //std::cout << setw(20) << tmpr;

                            //tmp = ab_p_values_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_p_derivatives_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_p_derivatives_in[k - 2];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_p_derivatives_in[k - 3];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_p_values_in[k];
                            //std::cout << setw(20) << tmp;
                            //std::cout << "\n";

                            ab_q_values_in[k] = ab_q_values_in[k - 1] + step * (23.0 * ab_q_derivatives_in[k - 1] - 16.0 * ab_q_derivatives_in[k - 2] + 5.0 * ab_q_derivatives_in[k - 3]) / 12.0;

                            if (ab_q_values_in[k] < 0.0) {
                                std::cout << "\n";
                                std::cout << setw(20) << "MASS BELOW ZERO";
                                std::cout << "\n";
                                break;
                            }

                            //std::cout << setw(20) << tmpr;

                            //tmp = ab_q_values_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_q_derivatives_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_q_derivatives_in[k - 2];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_q_derivatives_in[k - 3];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_q_values_in[k];
                            //std::cout << setw(20) << tmp;


                            ab_f_values_in[k] = ab_f_values_in[k - 1] + step * (23.0 * ab_f_derivatives_in[k - 1] - 16.0 * ab_f_derivatives_in[k - 2] + 5.0 * ab_f_derivatives_in[k - 3]) / 12.0;

                            if (ab_f_values_in[k] < 0.0) {
                                std::cout << "\n";
                                std::cout << setw(20) << "LUMINOSITY BELOW ZERO";
                                std::cout << setw(20) << tmpr;
                                tmp = ab_f_values_in[k - 1];
                                std::cout << setw(20) << tmp;
                                tmp = ab_f_derivatives_in[k - 1];
                                std::cout << setw(20) << tmp;
                                tmp = ab_f_derivatives_in[k - 2];
                                std::cout << setw(20) << tmp;
                                tmp = ab_f_derivatives_in[k - 3];
                                std::cout << setw(20) << tmp;
                                tmp = ab_f_values_in[k];
                                std::cout << setw(20) << tmp;
                                std::cout << "\n";
                                break;
                            }

                            //std::cout << setw(20) << tmpr;

                            //tmp = ab_f_values_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_f_derivatives_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_f_derivatives_in[k - 2];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_f_derivatives_in[k - 3];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_f_values_in[k];
                            //std::cout << setw(20) << tmp;


                            ab_t_values_in[k] = ab_t_values_in[k - 1] + step * (23.0 * ab_t_derivatives_in[k - 1] - 16.0 * ab_t_derivatives_in[k - 2] + 5.0 * ab_t_derivatives_in[k - 3]) / 12.0;

                            if (ab_t_values_in[k] < 0.0) {
                                std::cout << "\n";
                                std::cout << setw(20) << "TEMPERATURE BELOW ZERO";
                                std::cout << "\n";
                                break;
                            }

                            //std::cout << setw(20) << tmpr;

                            //tmp = ab_t_values_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_t_derivatives_in[k - 1];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_t_derivatives_in[k - 2];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_t_derivatives_in[k - 3];
                            //std::cout << setw(20) << tmp;
                            //tmp = ab_t_values_in[k];
                            //std::cout << setw(20) << tmp;


                            if (gas_equation_density == true) {
                                surface_density = (mu * ab_p_values_in[k] * proton_mass) / (boltzmann_constant * ab_t_values_in[k]);
                            }

                            if (gas_equation_density == false) {
                                // temporary estimate
                                surface_density = central_density / 60.0;
                            }

                            //tmp = surface_density;
                            //std::cout << setw(20) << tmp;
                            //std::cout << "\n";

                            // Schwarzschild page 81
                            epp_tentosix_over_t = 1000000.0 / ab_t_values_in[k];
                            epp_tentosix_over_t_to_twothirds = pow(epp_tentosix_over_t, 2.0 / 3.0);
                            epp_pressure_term = 2.5 * pow(10, 6) * surface_density * hydrogen_fraction * hydrogen_fraction;
                            epp_exponential_power = -33.8 * pow(epp_tentosix_over_t, 1.0 / 3.0);
                            epp_exponential_term = exp(epp_exponential_power);
                            surface_energy_generation = epp_pressure_term * epp_tentosix_over_t_to_twothirds * epp_exponential_term;
                            tmp = surface_energy_generation;
                            //std::cout << setw(20) << tmp;


                            ab_p_derivatives_in[k] = (surface_density * gravitational_constant * ab_q_values_in[k]) / pow(tmpr, 2);
                            ab_q_derivatives_in[k] = -surface_density * 4.0 * pi * pow(tmpr, 2);
                            ab_f_derivatives_in[k] = -surface_density * 4.0 * pi * pow(tmpr, 2) * surface_energy_generation;
                            ab_t_derivatives_in[k] = (3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (surface_opacity * surface_density / pow(ab_t_values_in[k], 3)) * (ab_f_values_in[k] / (4.0 * pi * pow(tmpr, 2)));
                        }

                        std::cout << "\n";
                        std::cout << "\n";
                        std::cout << "REVERSE ORDER RESULTS";
                        std::cout << "\n";
                        std::cout << "\n";

                        OutResults_File_mrl << "\n";
                        OutResults_File_mrl << "\n";
                        OutResults_File_mrl << "\n";
                        OutResults_File_mrl << setw(20) << "MASS = " << Model_Mass;
                        OutResults_File_mrl << setw(20) << "RADIUS = " << Model_Radius;
                        OutResults_File_mrl << setw(20) << "LUMINOSITY = " << Model_Luminosity;
                        OutResults_File_mrl << "\n";
                        OutResults_File_mrl << "\n";
                        OutResults_File_mrl << "\n";

                        for (int k = (ab_calculation_length - 10); k > numberofsteps; k--) {

                            if (gas_equation_density == true) {
                                surface_density = (mu * ab_p_values_in[k - 1] * proton_mass) / (boltzmann_constant * ab_t_values_in[k - 1]);
                            }

                            tmp = surface_density;
                            std::cout << setw(20) << tmp;
                            std::cout << "\n";

                            tmpr = Model_Radius - double(k) * step;
                            if (tmpr < 0.0) {
                                break;
                            }

                            std::cout << setw(20) << tmpr;
                            OutResults_File_mrl << setw(20) << tmpr;
                            OutResults_File_mrl << "P";

                            //if (k == ((ab_calculation_length / 2) - 1)) {
                            //    OutResults_File << Model_Radius;
                            //    OutResults_File << "          ";
                            //    OutResults_File << Model_Luminosity;
                            //    OutResults_File << "\n";
                            //    OutResults_File << ab_p_values_in[k];
                            //    OutResults_File << "\n";
                            //    OutResults_File << ab_q_values_in[k];
                            //    OutResults_File << "\n";
                            //    OutResults_File << ab_f_values_in[k];
                            //    OutResults_File << "\n";
                            //    OutResults_File << ab_t_values_in[k];
                            //    OutResults_File << "\n";
                            //    OutResults_File << "\n";
                            //}

                            tmp = ab_p_values_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "P";

                            tmp = ab_p_derivatives_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_p_derivatives_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_p_derivatives_in[k - 3];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_p_values_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "\n";

                            std::cout << "\n";

                            std::cout << setw(20) << tmpr;
                            OutResults_File_mrl << setw(20) << tmpr;
                            OutResults_File_mrl << "Q";

                            tmp = ab_q_values_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "Q";

                            tmp = ab_q_derivatives_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_q_derivatives_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_q_derivatives_in[k - 3];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_q_values_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "\n";
                            std::cout << "\n";

                            std::cout << setw(20) << tmpr;
                            OutResults_File_mrl << setw(20) << tmpr;
                            OutResults_File_mrl << "F";

                            tmp = ab_f_values_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "F";

                            tmp = ab_f_derivatives_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_f_derivatives_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_f_derivatives_in[k - 3];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_f_values_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "\n";
                            std::cout << "\n";

                            std::cout << setw(20) << tmpr;
                            OutResults_File_mrl << setw(20) << tmpr;
                            OutResults_File_mrl << "T";

                            tmp = ab_t_values_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "T";

                            tmp = ab_t_derivatives_in[k - 1];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_t_derivatives_in[k - 2];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;

                            tmp = ab_t_derivatives_in[k - 3];
                            std::cout << setw(20) << tmp;
                            OutResults_File_mrl << setw(20) << tmp;
                            tmp = ab_t_values_in[k - 2];
                            OutResults_File_mrl << setw(20) << tmp;
                            OutResults_File_mrl << "\n";

                            std::cout << setw(20) << tmp;
                            std::cout << "\n";

                        }
                    }
                }
            }
        }

        if (transformed) {
            //These are the variables of the first transformation - the unstarred variables
            //Schwarszchild page 116






        }
    }
}




//std::cout << setw(20) << "Transformed Pressure Variables Intermediate Values" << "    ";
//std::cout << setw(20) << C << "    ";
//std::cout << setw(20) << D << "    ";
//std::cout << setw(20) << p_zero << "    ";
//std::cout << setw(20) << pcentre_numerator << "    ";
//std::cout << setw(20) << pcentre_denominator << "    ";
//std::cout << setw(20) << pcentre << "    ";
//std::cout << setw(20) << pstarc << "    ";
//std::cout << "\n";
//std::cout << "\n";
//std::cout << "\n";

//std::cout << setw(20) << "Transformed Distance Variables Intermediate Values" << "    ";
//std::cout << setw(20) << x_zero << "    ";
//std::cout << "\n";
//std::cout << "\n";
//std::cout << "\n";


//if (tmpr > zero_pressure_distance) {
   //    std::cout << "\n";
   //    std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
   //    std::cout << "\n";
   //    break;
   //}

   // Now get the new first order derivatives



   //if (tmpr < Layer_Width)
   //{
   //    H_abund =  AH[0] + BH[0] * tmpr * pow(10, -10);
   //    He_abund = AHe[0] + BHe[0] * tmpr * pow(10, -10);
   //}

   //if (tmpr > Layer_Width && 2.0*Layer_Width > tmpr)
   //{
   //    H_abund = AH[1] + BH[1] * tmpr * pow(10, -10);
   //    He_abund = AHe[1] + BHe[1] * tmpr * pow(10, -10);
   //}

   //if (tmpr > 2.0*Layer_Width && 3.0*Layer_Width > tmpr)
   //{
   //    H_abund = AH[2] + BH[2] * tmpr * pow(10, -10);
   //    He_abund = AHe[2] + BHe[2] * tmpr * pow(10, -10);
   //}

   //if (tmpr > 3.0*Layer_Width)
   //{
   //    H_abund = AH[3] + BH[3] * tmpr * pow(10, -10);
   //    He_abund = AHe[3] + BHe[3] * tmpr * pow(10, -10);
   //}

   //Hv_abund = 1.0 - H_abund - He_abund;

   //mu_reciprocal = 2.0 * H_abund + 0.75 * He_abund + 0.5 * Hv_abund;
   //mu = 1.0 / mu_reciprocal;

   //std::cout << setw(20) << "ABUNDANCE VALUE";
   //std::cout << setw(20) << mu;
   //std::cout << "\n";

   //new_density = (mu * ab_p_values[k] * proton_mass) / (boltzmann_constant * ab_t_values[k]);
   ////Opacity and energy generation also need updates



