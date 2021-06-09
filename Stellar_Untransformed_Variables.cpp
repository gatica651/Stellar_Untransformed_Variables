




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
// Calculates a few starting values at the stellar centre using a constant density approximation
// Inputs these starting values into an Adams-Bashforth integration algorithm and calculates outwards towards the surface of the star
// Calculates a few starting values at the stellar surface using an approximation where mass and luminosity are assumed constant and at their full values
// Inputs these starting values into the same Adams-Bashforth integration algorithm and calculates inwards towards the centre of the star
// There is code for simple extrapolation below that is not used and could maybe be removed
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
// Do evolutionary calculations where the star develops in time and its chemical composition changes
// Apply the Henyey method of stellar structure calculation
// Do calculations with the stellar mass (instead of distance) as the independent variable.
//
// This code is not guaranteed to be bug-free nor to be in the most elegant form.
// Please use and develop the uploaded Papers 1-3 and this code as you wish but credit me please : 
//
// Peter Gardner : peter.gardner13@btopenworld.com, UK 07495 416625, also see my LinkedIn profile.
//
// The original uploaded papers and code are read-only.
// 





#include <iostream>
#include <iomanip>
using namespace std;

#define calculation_length 200

class extrapolated_values {
public:
    //static const int calculation_length = 200;
    double calculated_first[calculation_length];
    double calculated_second[calculation_length];
    double calculated_third[calculation_length];
    double calculated_fourth[calculation_length];
    double calculated_values[calculation_length];
};

class form_difference_table {
public:
    void create_difference_table(int numberofsteps, double start_values[], double first_differences[], double second_differences[], double third_differences[], double fourth_differences[])
    {

        double tmp = 0.0;

        for (int m = 1; m < numberofsteps; m++)
        {
            //std::cout << setw(20) << m;
            tmp = start_values[numberofsteps - m] - start_values[numberofsteps - 1 - m];
            first_differences[numberofsteps - 1 - m] = tmp;
            //std::cout << setw(20) << tmp;
        }

        std::cout << "\n";

        for (int m = 2; m < numberofsteps; m++)
        {
            //std::cout << setw(20) << m;
            tmp = first_differences[numberofsteps - m] - first_differences[numberofsteps - 1 - m];
            second_differences[numberofsteps - 1 - m] = tmp;
            //std::cout << setw(20) << tmp;
        }

        std::cout << "\n";

        for (int m = 3; m < numberofsteps; m++)
        {
            //std::cout << setw(20) << m;
            tmp = second_differences[numberofsteps - m] - second_differences[numberofsteps - 1 - m];
            third_differences[numberofsteps - 1 - m] = tmp;
            //std::cout << setw(20) << tmp;
        }

        std::cout << "\n";

        for (int m = 4; m < numberofsteps; m++)
        {
            //std::cout << setw(20) << m;
            tmp = third_differences[numberofsteps - m] - third_differences[numberofsteps - 1 - m];
            fourth_differences[numberofsteps - 1 - m] = tmp;
            //std::cout << setw(20) << tmp;
        }

    }


    void output_difference_table(int numberofsteps, double step, double start_values[], double first_differences[], double second_differences[], double third_differences[], double fourth_differences[])
    {

        double tmp = 0.0;
        double r = 0.0;

        for (int i = 1; i < 2 * numberofsteps; i++)
        {
            switch (i) {
            case 1:
                r = (double(i) - 1) * step;
                std::cout << setw(20) << r;
                tmp = start_values[0];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 2:
                std::cout << setw(20) << "    ";
                std::cout << setw(20) << "    ";
                tmp = first_differences[0];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 3:
                r = (double(i) - 2) * step;
                std::cout << setw(20) << r;
                tmp = start_values[1];
                std::cout << setw(20) << tmp;
                std::cout << setw(20) << "    ";
                tmp = second_differences[0];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 4:
                std::cout << setw(20) << "    ";
                std::cout << setw(20) << "    ";
                tmp = first_differences[1];
                std::cout << setw(20) << tmp;
                std::cout << setw(20) << "    ";
                tmp = third_differences[0];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 5:
                r = (double(i) - 3) * step;
                std::cout << setw(20) << r;
                tmp = start_values[2];
                std::cout << setw(20) << tmp;
                std::cout << setw(20) << "    ";
                tmp = second_differences[1];
                std::cout << setw(20) << tmp;
                std::cout << setw(20) << "    ";
                tmp = fourth_differences[0];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 6:
                std::cout << setw(20) << "    ";
                std::cout << setw(20) << "    ";
                tmp = first_differences[2];
                std::cout << setw(20) << tmp;
                std::cout << setw(20) << "    ";
                tmp = third_differences[1];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 7:
                r = (double(i) - 4) * step;
                std::cout << setw(20) << r;
                tmp = start_values[3];
                std::cout << setw(20) << tmp;
                std::cout << setw(20) << "    ";
                tmp = second_differences[2];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 8:
                std::cout << setw(20) << "    ";
                std::cout << setw(20) << "    ";
                tmp = first_differences[3];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            case 9:
                r = (double(i) - 5) * step;
                std::cout << setw(20) << r;
                tmp = start_values[4];
                std::cout << setw(20) << tmp;
                std::cout << "\n";
                break;
            }
        }


    }

};



int main(int argc, char* argv []) {


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

    double mu_reciprocal = 2.0 * hydrogen_fraction + 0.75 * helium_fraction + 0.5 * heavy_fraction;
    double mu = 1.0 / mu_reciprocal;                                                                      //Schwarszchild equation 8.2
    const double Central_Pressure = 6.0 * pow(10, 15);
    const double Central_Temperature = 1.56 * pow(10, 7);
    const double Total_Mass = 2.0 * pow(10, 33);                                                           //Wikipedia                           (grams)
    const double Radius = 7.0 * pow(10, 10);                                                               //Wikipedia                           (cms)
    //const double Total_Luminosity = 6.0 * pow(10, 33);   //Schwarzschild page 42
    const double Total_Luminosity = 3.0 * pow(10, 33);
    const double t_over_tentosix = 17.0;                                                                  //Schwarzschild page 83, Table 10.1
    const double epsilon1 = 0.00001;                                                                       //Schwarzschild page 83, Table 10.1
    const double nu = 4.0;                                                                                //Schwarzschild page 83, Table 10.1
    const double gbart = 1.0;                                                                             // Schwarzschild equation 9.16
    const double gffbar = 1.0;                                                                            // Schwarzschild equation 9.18
    bool gas_equation_density = false;
    double central_density = 0.0;
    double chi_zero = 0.0;
    double chi_bf = 0.0;
    double chi_ff = 0.0;
    double central_opacity = 0.0;

    if (gas_equation_density == true) {
        std::cout << "\n";
        central_density = (mu * Central_Pressure * proton_mass) / (boltzmann_constant * Central_Temperature);
        std::cout << "\n";
        std::cout << "\n";
    }


    if (gas_equation_density == false) {
        std::cout << "\n";
        central_density = 3.0;
        std::cout << "\n";
        std::cout << "\n";
    }

    chi_bf = (4.34 * pow(10, 25) * gbart * heavy_fraction * (1.0 + hydrogen_fraction));                              // Schwarzschild equation 9.16
    chi_ff = (3.68 * pow(10, 22) * gffbar * (hydrogen_fraction + helium_fraction) * (1.0 + hydrogen_fraction));      // Schwarzschild equation 9.18
    chi_zero = chi_bf + chi_ff;

    central_opacity = chi_zero*central_density/pow(Central_Temperature, 3.5);

    const double central_energy_generation = epsilon1 * central_density * hydrogen_fraction * hydrogen_fraction * pow(t_over_tentosix, nu);  //proton-proton only

    std::cout << "Total Mass =  " << Total_Mass;
    std::cout << "\n";
    std::cout << "Total Luminosity =  " << Total_Luminosity;
    std::cout << "\n";
    std::cout << "Composition =  " << mu;
    std::cout << "\n";
    std::cout << "Radius =  " << Radius;
    std::cout << "\n";
    std::cout << "Central Pressure =  " << Central_Pressure;
    std::cout << "\n";
    std::cout << "Central Temperature =  " << Central_Temperature;
    std::cout << "\n";
    std::cout << "Central Density =  " << central_density;
    std::cout << "\n";
    std::cout << "Central Energy Generation =  " << central_energy_generation;
    std::cout << "\n";
    std::cout << "Chi_bf =  " << chi_bf;
    std::cout << "\n";
    std::cout << "Chi_ff =  " << chi_ff;
    std::cout << "\n";
    std::cout << "Central Opacity =  " << central_opacity;
    std::cout << "\n";
    std::cout << "\n";

    //THE VARIABLES OF THE DIFFERENTIAL EQUATIONS

    double r = 0.0, from_surface;
    double rfactor, temperature_constants_factor;
    double pressure_fundamental_factors, pressure_stellar_factors, pressure_all_factors;
    double pressure, mass, luminosity, temperature;
    double zero_pressure_distance = 0.0;
    bool zero_pressure_reached = false;
    int zero_pressure_j = 0;

    // VARIABLES OF THE ALGORITHM

    const int numberofsteps = 5;
    bool extrapolation = false;
    bool adams_bashforth = true;
    bool untransformed = true;
    bool transformed = false;
    bool linear_density_decrease = true;
    double step = 0.0;
    double step_size = 0.02;
    if (untransformed) {
        step = 14.0 * pow(10, 8);
    } else {                                                //step size for untransformed variableselse {
        step = 0.1;                                        //step size for transformed variables
    }
    double new_density = 0.0;
    double density = 0.0;
    

    // VARIABLES OF THE ALGORITHM (OUTWARDS SOLUTION)

    double p_maclaurinterms[numberofsteps];
    double q_maclaurinterms[numberofsteps];
    double f_maclaurinterms[numberofsteps];
    double t_maclaurinterms[numberofsteps];

    double p_centre_derivatives[numberofsteps];
    double q_centre_derivatives[numberofsteps];
    double f_centre_derivatives[numberofsteps];
    double t_centre_derivatives[numberofsteps];

    double p_maclaurinterms_first[numberofsteps - 1];
    double p_maclaurinterms_second[numberofsteps - 2];
    double p_maclaurinterms_third[numberofsteps - 3];
    double p_maclaurinterms_fourth[numberofsteps - 4];

    double q_maclaurinterms_first[numberofsteps - 1];
    double q_maclaurinterms_second[numberofsteps - 2];
    double q_maclaurinterms_third[numberofsteps - 3];
    double q_maclaurinterms_fourth[numberofsteps - 4];

    double f_maclaurinterms_first[numberofsteps - 1];
    double f_maclaurinterms_second[numberofsteps - 2];
    double f_maclaurinterms_third[numberofsteps - 3];
    double f_maclaurinterms_fourth[numberofsteps - 4];

    double t_maclaurinterms_first[numberofsteps - 1];
    double t_maclaurinterms_second[numberofsteps - 2];
    double t_maclaurinterms_third[numberofsteps - 3];
    double t_maclaurinterms_fourth[numberofsteps - 4];

    // VARIABLES OF THE ALGORITHM (INWARDS SOLUTION)

    double p_surface[calculation_length];
    double q_surface[calculation_length];
    double f_surface[calculation_length];
    double t_surface[calculation_length];

    double p_surface_derivatives[numberofsteps];
    double q_surface_derivatives[numberofsteps];
    double f_surface_derivatives[numberofsteps];
    double t_surface_derivatives[numberofsteps];

    double p_surface_first[numberofsteps - 1];
    double p_surface_second[numberofsteps - 2];
    double p_surface_third[numberofsteps - 3];
    double p_surface_fourth[numberofsteps - 4];

    double q_surface_first[numberofsteps - 1];
    double q_surface_second[numberofsteps - 2];
    double q_surface_third[numberofsteps - 3];
    double q_surface_fourth[numberofsteps - 4];

    double f_surface_first[numberofsteps - 1];
    double f_surface_second[numberofsteps - 2];
    double f_surface_third[numberofsteps - 3];
    double f_surface_fourth[numberofsteps - 4];

    double t_surface_first[numberofsteps - 1];
    double t_surface_second[numberofsteps - 2];
    double t_surface_third[numberofsteps - 3];
    double t_surface_fourth[numberofsteps - 4];

    // VARIABLES OF THE ALGORITHM (TRANSFORMED VARIABLES)

    double xstar, pstar, qstar, fstar, tstar;  

    // Constant C : Schwarszchild equation 13.9
    double c1 = 3.0 / (4.0 * radiation_density_constant * speed_of_light);
    double c2 = boltzmann_constant / (proton_mass * gravitational_constant);
    double c3 = pow(c2, 7.5);
    double c4 = (1.0) / (pow(4.0 * pi, 3));
    double c5 = central_opacity / (pow(mu, 7.5));
    double c6 = Total_Luminosity * pow(Radius, 0.5) / pow(Total_Mass, 5.5);
    double C = c1 * c3 * c4 * c5 * c6;

    // Constant D : Schwarszchild equation 13.10
    double d1 = proton_mass * gravitational_constant / boltzmann_constant;
    double d2 = pow(d1, nu)/(4.0*pi);
    double d3 = central_energy_generation * pow(mu, nu);
    double d4 = pow(Total_Mass, nu + 2) / (Total_Luminosity * pow(Radius, nu + 3));
    double D = d2 * d3 * d4;

    // The starred variables : Schwarszchild equations 13.13 
    double tc = (boltzmann_constant * Radius * Central_Temperature) / (mu * proton_mass * gravitational_constant * Total_Mass);
    double p_zero_cubed = pow(tc, 9.5 - nu) / (C * D);
    //double p_zero = pow(p_zero_cubed, 1.0/3.0);
    double p_zero = cbrt(p_zero_cubed);
    double pcentre_numerator = 4.0 * pi * pow(Radius, 4);
    double pcentre_denominator =  gravitational_constant * pow(Total_Mass, 2);
    double pcentre = (pcentre_numerator / pcentre_denominator)*Central_Pressure;
    //The value of pstar at the centre
    double pstarc = pcentre / p_zero;

    double x_zero_1 = pow(Central_Temperature, 9.5 - nu + 2) / (C * D);
    double x_zero = pow(x_zero_1, 0.5) / pow(p_zero, 2);

    //Atomic abundance profile (H, He, heavies)
    const int NumberOfLayers = 4;
    const double Layer_Width = Radius / double(NumberOfLayers);
    double AH[NumberOfLayers] = {0.08, 0.53, 0.20, 0.06};       
    double BH[NumberOfLayers] = {0.34, 0.40, 0.80, 0.95};           // apply factor of pow(10, -10)
    double AHe[NumberOfLayers] = {0.08, 0.50, 0.20, 0.06 };       
    double BHe[NumberOfLayers] = {-0.64, -0.58, -0.20, -0.05 };     // apply factor of pow(10, -10)
    double H_abund = 0.0, He_abund = 0.0, Hv_abund = 0.0;

    double lhs = 0.0;
    double rhs = 0.0;

    std::cout << setw(20) << "Transformed Pressure Variables Intermediate Values" << "    ";
    std::cout << setw(20) << C << "    ";
    std::cout << setw(20) << D << "    ";
    std::cout << setw(20) << p_zero << "    ";
    std::cout << setw(20) << pcentre_numerator << "    ";
    std::cout << setw(20) << pcentre_denominator << "    ";
    std::cout << setw(20) << pcentre << "    ";
    std::cout << setw(20) << pstarc << "    ";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout << setw(20) << "Transformed Distance Variables Intermediate Values" << "    ";
    std::cout << setw(20) << x_zero << "    ";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    double tmp, tmpr;

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

    if (linear_density_decrease==true) {
        std::cout << "LINEAR DENSITY DECREASE";
    }
    if (linear_density_decrease==false) {
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

        for (int n = 0; n < numberofsteps; n++) {
            r = double(n) * step;
            density = central_density;

            if (linear_density_decrease==true) {
                density = central_density - central_density * double(n) * step_size;
            }

            // Schwarchild page 114

            pressure = Central_Pressure - (2.0 * pi * gravitational_constant * pow(density, 2) * pow(r, 2)) / 3.0;
            mass = (4.0 * pi * density * pow(r, 3)) / 3.0;
            luminosity = (4.0 * pi * density * central_energy_generation * pow(r, 3)) / 3.0;
            temperature = Central_Temperature - (central_opacity * central_energy_generation * pow(density, 2) * pow(r, 2)) / (8.0 * radiation_density_constant * speed_of_light * pow(Central_Temperature, 3));

            p_maclaurinterms[n] = pressure;
            q_maclaurinterms[n] = mass;
            f_maclaurinterms[n] = luminosity;
            t_maclaurinterms[n] = temperature;

            
            if (n == 0) {
                p_centre_derivatives[0] = 0.0;
                q_centre_derivatives[0] = 0.0;
                f_centre_derivatives[0] = 0.0;
                t_centre_derivatives[0] = 0.0;
            }

            if (n > 0) {
                p_centre_derivatives[n] = -(density * gravitational_constant * q_maclaurinterms[n]) / pow(r, 2);
                q_centre_derivatives[n] = density * 4.0 * pi * pow(r, 2);
                f_centre_derivatives[n] = density * 4.0 * pi * pow(r, 2) * central_energy_generation;
                t_centre_derivatives[n] = -(3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (central_opacity * density / pow(t_maclaurinterms[n], 3)) * (f_maclaurinterms[n] / (4.0 * pi * pow(r, 2)));
            }

            std::cout << setw(10) << r << "    ";
            std::cout << setw(10) << density << "    ";
            std::cout << setw(10) << pressure << "    ";
            std::cout << setw(10) << mass << "    ";
            std::cout << setw(10) << luminosity << "    ";
            std::cout << setw(10) << temperature << "    ";
            std::cout << setw(10) << p_centre_derivatives[n] << "    ";
            std::cout << setw(10) << q_centre_derivatives[n] << "    ";
            std::cout << setw(10) << f_centre_derivatives[n] << "    ";
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

            pstar = pstarc - pow(pstarc, 2) * pow(xstar, 2) / 6.0;
            qstar = pstarc * pow(xstar, 3) / 3.0;
            fstar = pow(pstarc, 2) * pow(xstar, 3) / 3.0;
            tstar = 1.0 - pow(pstarc, 4) * pow(xstar, 2) / 6.0;

            p_maclaurinterms[n] = pstar;
            q_maclaurinterms[n] = qstar;
            f_maclaurinterms[n] = fstar;
            t_maclaurinterms[n] = tstar;

            std::cout << setw(20) << xstar << "    ";
            std::cout << setw(20) << pstar << "    ";
            std::cout << setw(20) << qstar << "    ";
            std::cout << setw(20) << fstar << "    ";
            std::cout << setw(20) << tstar << "    ";
            std::cout << "\n";


        }
    }
        
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    // Extrapolation Integration

    std::cout << "CREATING DIFFERENCE TABLES";

    std::cout << "\n";
    std::cout << "\n";

    std::cout << "PRESSURE";

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";


    form_difference_table difference_table;

    difference_table.create_difference_table(numberofsteps, p_maclaurinterms, p_maclaurinterms_first, p_maclaurinterms_second, p_maclaurinterms_third, p_maclaurinterms_fourth);
    difference_table.output_difference_table(numberofsteps, step, p_maclaurinterms, p_maclaurinterms_first, p_maclaurinterms_second, p_maclaurinterms_third, p_maclaurinterms_fourth);
    std::cout << "\n";
    std::cout << "\n";

    std::cout << "MASS";

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    difference_table.create_difference_table(numberofsteps, q_maclaurinterms, q_maclaurinterms_first, q_maclaurinterms_second, q_maclaurinterms_third, q_maclaurinterms_fourth);
    difference_table.output_difference_table(numberofsteps, step, q_maclaurinterms, q_maclaurinterms_first, q_maclaurinterms_second, q_maclaurinterms_third, q_maclaurinterms_fourth);
    std::cout << "\n";
    std::cout << "\n";

    std::cout << "LUMINOSITY";

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    difference_table.create_difference_table(numberofsteps, f_maclaurinterms, f_maclaurinterms_first, f_maclaurinterms_second, f_maclaurinterms_third, f_maclaurinterms_fourth);
    difference_table.output_difference_table(numberofsteps, step, f_maclaurinterms, f_maclaurinterms_first, f_maclaurinterms_second, f_maclaurinterms_third, f_maclaurinterms_fourth);
    std::cout << "\n";
    std::cout << "\n";

    std::cout << "TEMPERATURE";

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    difference_table.create_difference_table(numberofsteps, t_maclaurinterms, t_maclaurinterms_first, t_maclaurinterms_second, t_maclaurinterms_third, t_maclaurinterms_fourth);
    difference_table.output_difference_table(numberofsteps, step, t_maclaurinterms, t_maclaurinterms_first, t_maclaurinterms_second, t_maclaurinterms_third, t_maclaurinterms_fourth);

    std::cout << "\n";
    std::cout << "\n";

    if (extrapolation) {
        std::cout << "INTEGRATING BY EXTRAPOLATION";

        std::cout << "\n";
        std::cout << "\n";
        std::cout << "PRESSURE";
        std::cout << "\n";
        std::cout << "\n";

        extrapolated_values extrapolated_pressure_values;

        extrapolated_pressure_values.calculated_fourth[0] = p_maclaurinterms_fourth[0];
        extrapolated_pressure_values.calculated_third[0] = p_maclaurinterms_third[1] + extrapolated_pressure_values.calculated_fourth[0];
        extrapolated_pressure_values.calculated_second[0] = p_maclaurinterms_second[2] + extrapolated_pressure_values.calculated_third[0];
        extrapolated_pressure_values.calculated_first[0] = p_maclaurinterms_first[3] + extrapolated_pressure_values.calculated_second[0];
        extrapolated_pressure_values.calculated_values[0] = p_maclaurinterms[numberofsteps - 1] + extrapolated_pressure_values.calculated_first[0];

        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_pressure_values.calculated_fourth[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_pressure_values.calculated_third[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_pressure_values.calculated_second[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_pressure_values.calculated_first[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        r = (numberofsteps)*step;
        std::cout << setw(20) << r;
        tmp = extrapolated_pressure_values.calculated_values[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";

        for (int j = 1; j < calculation_length; j++) {
            extrapolated_pressure_values.calculated_fourth[j] = extrapolated_pressure_values.calculated_fourth[j - 1];
            extrapolated_pressure_values.calculated_third[j] = extrapolated_pressure_values.calculated_third[j - 1] + extrapolated_pressure_values.calculated_fourth[j - 1];
            extrapolated_pressure_values.calculated_second[j] = extrapolated_pressure_values.calculated_second[j - 1] + extrapolated_pressure_values.calculated_third[j - 1];
            extrapolated_pressure_values.calculated_first[j] = extrapolated_pressure_values.calculated_first[j - 1] + extrapolated_pressure_values.calculated_second[j - 1];
            extrapolated_pressure_values.calculated_values[j] = extrapolated_pressure_values.calculated_values[j - 1] + extrapolated_pressure_values.calculated_first[j - 1];
            if (extrapolated_pressure_values.calculated_values[j] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "PRESSURE BELOW ZERO";
                std::cout << "\n";
                zero_pressure_reached = true;
                zero_pressure_distance = double(j + 5) * step;
                zero_pressure_j = j;
                break;
            }

            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_pressure_values.calculated_fourth[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_pressure_values.calculated_third[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_pressure_values.calculated_second[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_pressure_values.calculated_first[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            r = (numberofsteps + j) * step;
            std::cout << setw(20) << r;
            tmp = extrapolated_pressure_values.calculated_values[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
        }

        std::cout << "MASS";

        std::cout << "\n";
        std::cout << "\n";

        extrapolated_values extrapolated_mass_values;

        extrapolated_mass_values.calculated_fourth[0] = q_maclaurinterms_fourth[0];
        extrapolated_mass_values.calculated_third[0] = q_maclaurinterms_third[1] + extrapolated_mass_values.calculated_fourth[0];
        extrapolated_mass_values.calculated_second[0] = q_maclaurinterms_second[2] + extrapolated_mass_values.calculated_third[0];
        extrapolated_mass_values.calculated_first[0] = q_maclaurinterms_first[3] + extrapolated_mass_values.calculated_second[0];
        extrapolated_mass_values.calculated_values[0] = q_maclaurinterms[numberofsteps - 1] + extrapolated_mass_values.calculated_first[0];

        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_mass_values.calculated_fourth[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_mass_values.calculated_third[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_mass_values.calculated_second[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_mass_values.calculated_first[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        r = (numberofsteps)*step;
        std::cout << setw(20) << r;
        tmp = extrapolated_mass_values.calculated_values[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";

        for (int j = 1; j < calculation_length; j++) {
            extrapolated_mass_values.calculated_fourth[j] = extrapolated_mass_values.calculated_fourth[j - 1];
            extrapolated_mass_values.calculated_third[j] = extrapolated_mass_values.calculated_third[j - 1] + extrapolated_mass_values.calculated_fourth[j - 1];
            extrapolated_mass_values.calculated_second[j] = extrapolated_mass_values.calculated_second[j - 1] + extrapolated_mass_values.calculated_third[j - 1];
            extrapolated_mass_values.calculated_first[j] = extrapolated_mass_values.calculated_first[j - 1] + extrapolated_mass_values.calculated_second[j - 1];
            extrapolated_mass_values.calculated_values[j] = extrapolated_mass_values.calculated_values[j - 1] + extrapolated_mass_values.calculated_first[j - 1];
            if (extrapolated_mass_values.calculated_values[j] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "MASS BELOW ZERO";
                std::cout << "\n";
                break;
            }


            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_mass_values.calculated_fourth[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_mass_values.calculated_third[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_mass_values.calculated_second[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_mass_values.calculated_first[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            r = (numberofsteps + j) * step;
            std::cout << setw(20) << r;
            tmp = extrapolated_mass_values.calculated_values[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";

            if (r > zero_pressure_distance) {
                std::cout << "\n";
                std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
                std::cout << "\n";
                break;
            }
        }


         std::cout << "LUMINOSITY";

         std::cout << "\n";
         std::cout << "\n";

         extrapolated_values extrapolated_luminosity_values;

         extrapolated_luminosity_values.calculated_fourth[0] = f_maclaurinterms_fourth[0];
         extrapolated_luminosity_values.calculated_third[0] = f_maclaurinterms_third[1] + extrapolated_luminosity_values.calculated_fourth[0];
         extrapolated_luminosity_values.calculated_second[0] = f_maclaurinterms_second[2] + extrapolated_luminosity_values.calculated_third[0];
         extrapolated_luminosity_values.calculated_first[0] = f_maclaurinterms_first[3] + extrapolated_luminosity_values.calculated_second[0];
         extrapolated_luminosity_values.calculated_values[0] = f_maclaurinterms[numberofsteps - 1] + extrapolated_luminosity_values.calculated_first[0];

        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_luminosity_values.calculated_fourth[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_luminosity_values.calculated_third[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_luminosity_values.calculated_second[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_luminosity_values.calculated_first[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        r = (numberofsteps)*step;
        std::cout << setw(20) << r;
        tmp = extrapolated_luminosity_values.calculated_values[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";

        for (int j = 1; j < calculation_length; j++) {
            extrapolated_luminosity_values.calculated_fourth[j] = extrapolated_luminosity_values.calculated_fourth[j - 1];
            extrapolated_luminosity_values.calculated_third[j] = extrapolated_luminosity_values.calculated_third[j - 1] + extrapolated_luminosity_values.calculated_fourth[j - 1];
            extrapolated_luminosity_values.calculated_second[j] = extrapolated_luminosity_values.calculated_second[j - 1] + extrapolated_luminosity_values.calculated_third[j - 1];
            extrapolated_luminosity_values.calculated_first[j] = extrapolated_luminosity_values.calculated_first[j - 1] + extrapolated_luminosity_values.calculated_second[j - 1];
            extrapolated_luminosity_values.calculated_values[j] = extrapolated_luminosity_values.calculated_values[j - 1] + extrapolated_luminosity_values.calculated_first[j - 1];
            if (extrapolated_luminosity_values.calculated_values[j] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "LUMINOSITY BELOW ZERO";
                std::cout << "\n";
                break;
            }


            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_luminosity_values.calculated_fourth[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_luminosity_values.calculated_third[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_luminosity_values.calculated_second[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_luminosity_values.calculated_first[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            r = (numberofsteps + j) * step;
            std::cout << setw(20) << r;
            tmp = extrapolated_luminosity_values.calculated_values[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";

            if (r > zero_pressure_distance) {
                std::cout << "\n";
                std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
                std::cout << "\n";
                break;
            }
        }

        std::cout << "TEMPERATURE";

        std::cout << "\n";
        std::cout << "\n";

        extrapolated_values extrapolated_temperature_values;

        extrapolated_temperature_values.calculated_fourth[0] = t_maclaurinterms_fourth[0];
        extrapolated_temperature_values.calculated_third[0] = t_maclaurinterms_third[1] + extrapolated_temperature_values.calculated_fourth[0];
        extrapolated_temperature_values.calculated_second[0] = t_maclaurinterms_second[2] + extrapolated_temperature_values.calculated_third[0];
        extrapolated_temperature_values.calculated_first[0] = t_maclaurinterms_first[3] + extrapolated_temperature_values.calculated_second[0];
        extrapolated_temperature_values.calculated_values[0] = t_maclaurinterms[numberofsteps - 1] + extrapolated_temperature_values.calculated_first[0];

        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_temperature_values.calculated_fourth[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //std::cout << setw(20) << "   ";
        //tmp = extrapolated_temperature_values.calculated_third[0];
        //std::cout << setw(20) << tmp;
        //std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_temperature_values.calculated_second[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        std::cout << setw(20) << "   ";
        std::cout << setw(20) << "   ";
        tmp = extrapolated_temperature_values.calculated_first[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";
        r = (numberofsteps)*step;
        std::cout << setw(20) << r;
        tmp = extrapolated_temperature_values.calculated_values[0];
        std::cout << setw(20) << tmp;
        std::cout << "\n";

        for (int j = 1; j < calculation_length; j++) {
            extrapolated_temperature_values.calculated_fourth[j] = extrapolated_temperature_values.calculated_fourth[j - 1];
            extrapolated_temperature_values.calculated_third[j] = extrapolated_temperature_values.calculated_third[j - 1] + extrapolated_temperature_values.calculated_fourth[j - 1];
            extrapolated_temperature_values.calculated_second[j] = extrapolated_temperature_values.calculated_second[j - 1] + extrapolated_temperature_values.calculated_third[j - 1];
            extrapolated_temperature_values.calculated_first[j] = extrapolated_temperature_values.calculated_first[j - 1] + extrapolated_temperature_values.calculated_second[j - 1];
            extrapolated_temperature_values.calculated_values[j] = extrapolated_temperature_values.calculated_values[j - 1] + extrapolated_temperature_values.calculated_first[j - 1];
            if (extrapolated_temperature_values.calculated_values[j] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "TEMPERATURE BELOW ZERO";
                std::cout << "\n";
                break;
            }


            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_temperature_values.calculated_fourth[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //std::cout << setw(20) << "   ";
            //tmp = extrapolated_temperature_values.calculated_third[j];
            //std::cout << setw(20) << tmp;
            //std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_temperature_values.calculated_second[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            std::cout << setw(20) << "   ";
            std::cout << setw(20) << "   ";
            tmp = extrapolated_temperature_values.calculated_first[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";
            r = (numberofsteps + j) * step;
            std::cout << setw(20) << r;
            tmp = extrapolated_temperature_values.calculated_values[j];
            std::cout << setw(20) << tmp;
            std::cout << "\n";

            if (r > zero_pressure_distance) {
                std::cout << "\n";
                std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
                std::cout << "\n";
                break;
            }
        }


        for (int j = 1; j < calculation_length; j++) {
            lhs = (2.0 * gravitational_constant * extrapolated_mass_values.calculated_values[j]) / (5.0 * extrapolated_pressure_values.calculated_values[j]);
            rhs = c1 * central_opacity * extrapolated_luminosity_values.calculated_values[j] / pow(extrapolated_temperature_values.calculated_values[j], 4);
            std::cout << "\n";
            std::cout << setw(20) << lhs;
            std::cout << setw(20) << rhs;
            std::cout << "\n";
            if (j == (zero_pressure_j-1)) {
                break;
            }
        }
    }

    if (adams_bashforth) {
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "INTEGRATING OUTWARD USING ADAMS-BASHFORTH METHOD";
        std::cout << "\n";

        const int ab_calculation_length = calculation_length + 5;   
        double ab_p_values[ab_calculation_length], ab_p_derivatives[ab_calculation_length];

        for (int k = 0; k < numberofsteps; k++) {
            ab_p_values[k] = p_maclaurinterms[k];
        }

        double ab_q_values[ab_calculation_length], ab_q_derivatives[ab_calculation_length];

        for (int k = 0; k < numberofsteps; k++) {
            ab_q_values[k] = q_maclaurinterms[k];
        }

        double ab_f_values[ab_calculation_length], ab_f_derivatives[ab_calculation_length];

        for (int k = 0; k < numberofsteps; k++) {
            ab_f_values[k] = f_maclaurinterms[k];
        }

        double ab_t_values[ab_calculation_length], ab_t_derivatives[ab_calculation_length];

        for (int k = 0; k < numberofsteps; k++) {
            ab_t_values[k] = t_maclaurinterms[k];
        }


        ab_p_derivatives[0] = 0.0;                  //value not needed
        ab_p_derivatives[1] = 0.0;                  //value not needed
        ab_p_derivatives[2] = p_centre_derivatives[2];
        ab_p_derivatives[3] = p_centre_derivatives[3];
        ab_p_derivatives[4] = p_centre_derivatives[4];
        ab_q_derivatives[0] = 0.0;                  //value not needed
        ab_q_derivatives[1] = 0.0;                  //value not needed
        ab_q_derivatives[2] = q_centre_derivatives[2];
        ab_q_derivatives[3] = q_centre_derivatives[3];
        ab_q_derivatives[4] = q_centre_derivatives[4];
        ab_f_derivatives[0] = 0.0;                  //value not needed
        ab_f_derivatives[1] = 0.0;                  //value not needed
        ab_f_derivatives[2] = f_centre_derivatives[2];
        ab_f_derivatives[3] = f_centre_derivatives[3];
        ab_f_derivatives[4] = f_centre_derivatives[4];
        ab_t_derivatives[0] = 0.0;                  //value not needed
        ab_t_derivatives[1] = 0.0;                  //value not needed
        ab_t_derivatives[2] = t_centre_derivatives[2];
        ab_t_derivatives[3] = t_centre_derivatives[3];
        ab_t_derivatives[4] = t_centre_derivatives[4];


        for (int k = 5; k < ab_calculation_length; k++) {


            //Adams-Bashforth third order predictor method

            tmpr = double(k) * step;
            std::cout << setw(20) << tmpr;

            ab_p_values[k] = ab_p_values[k - 1] + step * (23.0 * ab_p_derivatives[k - 1] - 16.0 * ab_p_derivatives[k - 2] + 5.0 * ab_p_derivatives[k - 3]) / 12.0;
            if (ab_p_values[k] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "PRESSURE BELOW ZERO";
                std::cout << "\n";
                break;
            }


            tmp = ab_p_values[k-1];
            std::cout << setw(20) << tmp;
            tmp = ab_p_derivatives[k - 1];
            std::cout << setw(20) << tmp;
            tmp = ab_p_derivatives[k - 2];
            std::cout << setw(20) << tmp;
            tmp = ab_p_derivatives[k - 3];
            std::cout << setw(20) << tmp;
            tmp = ab_p_values[k];
            std::cout << setw(20) << tmp;
            std::cout << "\n";

            ab_q_values[k] = ab_q_values[k - 1] + step * (23.0 * ab_q_derivatives[k - 1] - 16.0 * ab_q_derivatives[k - 2] + 5.0 * ab_q_derivatives[k - 3]) / 12.0;

            if (ab_q_values[k] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "MASS BELOW ZERO";
                    std::cout << "\n";
                break;
            }

            std::cout << setw(20) << tmpr;

            tmp = ab_q_values[k-1];
            std::cout << setw(20) << tmp;
            tmp = ab_q_derivatives[k - 1];
            std::cout << setw(20) << tmp;
            tmp = ab_q_derivatives[k - 2];
            std::cout << setw(20) << tmp;
            tmp = ab_q_derivatives[k - 3];
            std::cout << setw(20) << tmp;
            tmp = ab_q_values[k];
            std::cout << setw(20) << tmp;
            //if (tmpr > zero_pressure_distance) {
            //    std::cout << "\n";
            //    std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
            //    std::cout << "\n";
            //    break;
            //}
       
            ab_f_values[k] = ab_f_values[k - 1] + step * (23.0 * ab_f_derivatives[k - 1] - 16.0 * ab_f_derivatives[k - 2] + 5.0 * ab_f_derivatives[k - 3]) / 12.0;

            if (ab_f_values[k] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "LUMINOSITY BELOW ZERO";
                std::cout << "\n";
                break;
            }

            std::cout << setw(20) << tmpr;

            tmp = ab_f_values[k-1];
            std::cout << setw(20) << tmp;
            tmp = ab_f_derivatives[k - 1];
            std::cout << setw(20) << tmp;
            tmp = ab_f_derivatives[k - 2];
            std::cout << setw(20) << tmp;
            tmp = ab_f_derivatives[k - 3];
            std::cout << setw(20) << tmp;
            tmp = ab_f_values[k];
            std::cout << setw(20) << tmp;
            //if (tmpr > zero_pressure_distance) {
            //    std::cout << "\n";
            //    std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
            //    std::cout << "\n";
            //    break;
            //}
    
            ab_t_values[k] = ab_t_values[k - 1] + step * (23.0 * ab_t_derivatives[k - 1] - 16.0 * ab_t_derivatives[k - 2] + 5.0 * ab_t_derivatives[k - 3]) / 12.0;

            if (ab_t_values[k] < 0.0) {
                std::cout << "\n";
                std::cout << setw(20) << "TEMPERATURE BELOW ZERO";
                std::cout << "\n";
                break;
            }

            std::cout << setw(20) << tmpr;

            tmp = ab_t_values[k-1];
            std::cout << setw(20) << tmp;
            tmp = ab_t_derivatives[k - 1];
            std::cout << setw(20) << tmp;
            tmp = ab_t_derivatives[k - 2];
            std::cout << setw(20) << tmp;
            tmp = ab_t_derivatives[k - 3];
            std::cout << setw(20) << tmp;
            tmp = ab_t_values[k];
            std::cout << setw(20) << tmp;
   
            ab_p_derivatives[k] = - (density * gravitational_constant * ab_q_values[k]) / pow(tmpr, 2);
            ab_q_derivatives[k] = density * 4.0 * pi * pow(tmpr, 2);
            ab_f_derivatives[k] = density * 4.0 * pi * pow(tmpr, 2)*central_energy_generation;
            ab_t_derivatives[k] = -(3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (central_opacity * density / pow(ab_t_values[k], 3)) * (ab_f_values[k] / (4.0 * pi * pow(tmpr, 2)));
            std::cout << "\n";
        }
    }



    //Surface Values

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout << "SURFACE VALUES GOING IN";

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    if (untransformed) {

        std::cout << "UNTRANSFORMED VARIABLES";

        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";

        p_surface[0] = 0.0;
        q_surface[0] = Total_Mass;
        f_surface[0] = Total_Luminosity;
        t_surface[0] = 0.0;

        std::cout << setw(20) << Radius << "    ";
        std::cout << setw(20) << 0.0 << "    ";
        std::cout << setw(20) << Total_Mass << "    ";
        std::cout << setw(20) << Total_Luminosity << "    ";
        std::cout << setw(20) << 0.0 << "    ";
        std::cout << "\n";

        temperature_constants_factor = (mu * proton_mass * gravitational_constant * Total_Mass) / (4.25 * boltzmann_constant);        //radiative
        //constants_factor = (mu * proton_mass * gravitational_constant * mass) / (2.5 * boltzmann_constant);       //convective
        pressure_fundamental_factors = (2.0 / 8.5) * (4.0 * radiation_density_constant * speed_of_light / 3.0) * (boltzmann_constant / proton_mass);
        pressure_stellar_factors = (4.0 * pi * gravitational_constant) / (mu * chi_zero);
        pressure_all_factors = pow((pressure_fundamental_factors * pressure_stellar_factors) * (Total_Mass / Total_Luminosity), 0.5);


        std::cout << setw(20) << mu << "    ";
        std::cout << setw(20) << chi_zero << "    ";
        std::cout << setw(20) << temperature_constants_factor << "    ";
        std::cout << setw(20) << pressure_fundamental_factors << "    ";
        std::cout << setw(20) << pressure_stellar_factors << "    ";
        std::cout << setw(20) << pressure_all_factors << "    ";
        std::cout << "\n";
        std::cout << "\n";

        double surface_density = 0.0;
        double surface_pressure = 0.0;
        double surface_temperature = 0.0;
        double surface_opacity = 0.0;
        const double surface_t_over_tentosix = 5.0;                                                                      //Schwarzschild page 83, Table 10.1
        const double surface_epsilon1 = 0.0000001;                                                                       //Schwarzschild page 83, Table 10.1
        const double surface_nu = 6.0;                                                                                   //Schwarzschild page 83, Table 10.1
        double surface_energy_generation = 0.0;

        for (int n = 0; n < numberofsteps; n++) {
            from_surface = double(n) * step;
            r = Radius - from_surface;
            if (r < 0) {
                break;
            }
            rfactor = (1.0 / r) - (1.0 / Radius);

            surface_temperature = temperature_constants_factor * rfactor;
            surface_pressure = pressure_all_factors * pow(surface_temperature, 4.25);                 //radiative

            p_surface[n] = surface_pressure;
            q_surface[n] = Total_Mass;
            f_surface[n] = Total_Luminosity;
            t_surface[n] = surface_temperature;

            if (gas_equation_density == true) {
                surface_density = (mu * surface_pressure * proton_mass) / (boltzmann_constant * surface_temperature);
            }


            if (gas_equation_density == false) {
                // temporary estimate
                surface_density = central_density/10.0;
            }

            surface_opacity = chi_zero * surface_density / pow(surface_temperature, 3.5);
            surface_energy_generation = surface_epsilon1 * surface_density * hydrogen_fraction * hydrogen_fraction * pow(surface_t_over_tentosix, surface_nu);  //proton-proton only

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


            std::cout << setw(20) << r << "    ";
            std::cout << setw(20) << surface_density << "    ";
            std::cout << setw(20) << surface_pressure << "    ";
            std::cout << setw(20) << Total_Mass << "    ";
            std::cout << setw(20) << Total_Luminosity << "    ";
            std::cout << setw(20) << surface_temperature << "    ";
            std::cout << setw(10) << p_surface_derivatives[n] << "    ";
            std::cout << setw(10) << q_surface_derivatives[n] << "    ";
            std::cout << setw(10) << f_surface_derivatives[n] << "    ";
            std::cout << setw(10) << t_surface_derivatives[n] << "    ";
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

            const int ab_calculation_length = calculation_length + 5;
            double ab_p_values[ab_calculation_length], ab_p_derivatives[ab_calculation_length];

            for (int k = 0; k < numberofsteps; k++) {
                ab_p_values[k] = p_surface[k];
            }

            double ab_q_values[ab_calculation_length], ab_q_derivatives[ab_calculation_length];

            for (int k = 0; k < numberofsteps; k++) {
                ab_q_values[k] = q_surface[k];
            }

            double ab_f_values[ab_calculation_length], ab_f_derivatives[ab_calculation_length];

            for (int k = 0; k < numberofsteps; k++) {
                ab_f_values[k] = f_surface[k];
            }

            double ab_t_values[ab_calculation_length], ab_t_derivatives[ab_calculation_length];

            for (int k = 0; k < numberofsteps; k++) {
                ab_t_values[k] = t_surface[k];
            }

            ab_p_derivatives[0] = 0.0;                  //value not needed
            ab_p_derivatives[1] = 0.0;                  //value not needed
            ab_p_derivatives[2] = p_surface_derivatives[2];
            ab_p_derivatives[3] = p_surface_derivatives[3];
            ab_p_derivatives[4] = p_surface_derivatives[4];
            ab_q_derivatives[0] = 0.0;                  //value not needed
            ab_q_derivatives[1] = 0.0;                  //value not needed
            ab_q_derivatives[2] = q_surface_derivatives[2];
            ab_q_derivatives[3] = q_surface_derivatives[3];
            ab_q_derivatives[4] = q_surface_derivatives[4];
            ab_f_derivatives[0] = 0.0;                  //value not needed
            ab_f_derivatives[1] = 0.0;                  //value not needed
            ab_f_derivatives[2] = f_surface_derivatives[2];
            ab_f_derivatives[3] = f_surface_derivatives[3];
            ab_f_derivatives[4] = f_surface_derivatives[4];
            ab_t_derivatives[0] = 0.0;                  //value not needed
            ab_t_derivatives[1] = 0.0;                  //value not needed
            ab_t_derivatives[2] = t_surface_derivatives[2];
            ab_t_derivatives[3] = t_surface_derivatives[3];
            ab_t_derivatives[4] = t_surface_derivatives[4];

            for (int k = 5; k < ab_calculation_length; k++) {


                //Adams-Bashforth third order predictor method

                tmpr = Radius - double(k) * step;
                if (tmpr < 0.0) {
                    break;
                }
                std::cout << setw(20) << tmpr;

                ab_p_values[k] = ab_p_values[k - 1] + step * (23.0 * ab_p_derivatives[k - 1] - 16.0 * ab_p_derivatives[k - 2] + 5.0 * ab_p_derivatives[k - 3]) / 12.0;

                if (ab_p_values[k] < 0.0) {
                    std::cout << "\n";
                    std::cout << setw(20) << "PRESSURE BELOW ZERO";
                    std::cout << "\n";
                    break;
                }


                tmp = ab_p_values[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_p_derivatives[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_p_derivatives[k - 2];
                std::cout << setw(20) << tmp;
                tmp = ab_p_derivatives[k - 3];
                std::cout << setw(20) << tmp;
                tmp = ab_p_values[k];
                std::cout << setw(20) << tmp;
                std::cout << "\n";

                ab_q_values[k] = ab_q_values[k - 1] + step * (23.0 * ab_q_derivatives[k - 1] - 16.0 * ab_q_derivatives[k - 2] + 5.0 * ab_q_derivatives[k - 3]) / 12.0;

                if (ab_q_values[k] < 0.0) {
                    std::cout << "\n";
                    std::cout << setw(20) << "MASS BELOW ZERO";
                    std::cout << "\n";
                    break;
                }

                std::cout << setw(20) << tmpr;

                tmp = ab_q_values[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_q_derivatives[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_q_derivatives[k - 2];
                std::cout << setw(20) << tmp;
                tmp = ab_q_derivatives[k - 3];
                std::cout << setw(20) << tmp;
                tmp = ab_q_values[k];
                std::cout << setw(20) << tmp;
                //if (tmpr > zero_pressure_distance) {
                //    std::cout << "\n";
                //    std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
                //    std::cout << "\n";
                //    break;
                //}

                ab_f_values[k] = ab_f_values[k - 1] + step * (23.0 * ab_f_derivatives[k - 1] - 16.0 * ab_f_derivatives[k - 2] + 5.0 * ab_f_derivatives[k - 3]) / 12.0;

                if (ab_f_values[k] < 0.0) {
                    std::cout << "\n";
                    std::cout << setw(20) << "LUMINOSITY BELOW ZERO";
                    std::cout << "\n";
                    break;
                }

                std::cout << setw(20) << tmpr;

                tmp = ab_f_values[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_f_derivatives[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_f_derivatives[k - 2];
                std::cout << setw(20) << tmp;
                tmp = ab_f_derivatives[k - 3];
                std::cout << setw(20) << tmp;
                tmp = ab_f_values[k];
                std::cout << setw(20) << tmp;
                //if (tmpr > zero_pressure_distance) {
                //    std::cout << "\n";
                //    std::cout << setw(20) << "GONE BEYOND ZERO PRESSURE DISTANCE";
                //    std::cout << "\n";
                //    break;
                //}

                ab_t_values[k] = ab_t_values[k - 1] + step * (23.0 * ab_t_derivatives[k - 1] - 16.0 * ab_t_derivatives[k - 2] + 5.0 * ab_t_derivatives[k - 3]) / 12.0;

                if (ab_t_values[k] < 0.0) {
                    std::cout << "\n";
                    std::cout << setw(20) << "TEMPERATURE BELOW ZERO";
                    std::cout << "\n";
                    break;
                }

                std::cout << setw(20) << tmpr;

                tmp = ab_t_values[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_t_derivatives[k - 1];
                std::cout << setw(20) << tmp;
                tmp = ab_t_derivatives[k - 2];
                std::cout << setw(20) << tmp;
                tmp = ab_t_derivatives[k - 3];
                std::cout << setw(20) << tmp;
                tmp = ab_t_values[k];
                std::cout << setw(20) << tmp;

                ab_p_derivatives[k] = (surface_density * gravitational_constant * ab_q_values[k]) / pow(tmpr, 2);
                ab_q_derivatives[k] = -surface_density * 4.0 * pi * pow(tmpr, 2);
                ab_f_derivatives[k] = -surface_density * 4.0 * pi * pow(tmpr, 2) * surface_energy_generation;
                ab_t_derivatives[k] = (3.0 / (4.0 * radiation_density_constant * speed_of_light)) * (surface_opacity * surface_density / pow(ab_t_values[k], 3)) * (ab_f_values[k] / (4.0 * pi * pow(tmpr, 2)));
                std::cout << "\n";
            }
        }


    }

    if (transformed) {
        //These are the variables of the first transformation - the unstarred variables
        //Schwarszchild page 116






    }

}




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

