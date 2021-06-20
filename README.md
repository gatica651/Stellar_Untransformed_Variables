# Stellar_Untransformed_Variables
                               C++ CODE FOR THE NUMERICAL CALCULATION OF STELLAR INTERIORS
                               
This project contains C++ code to calculate the structure of main sequence stars.

The code is derived from equations in Papers 2 and 3 of this series :

+ Paper 1 : TechnicalNote_StellarStructure         (this is TBD, not uploaded yet)
+ Paper 2 : TechnicalNote_NumericalIntegration     (this is uploaded and substantially complete)
+ Paper 3 : TechnicalNote_CplusplusCode            (this is uploaded but lacks sections on transformed variables and Henyey calculations)

Papers 2 and 3 are also in the project as can be seen in the listing above.

The code does essentially four things :

+ Calculates a few starting values at the stellar centre using a modified constant density approximation
+ Inputs these starting values into an Adams-Bashforth integration algorithm and calculates outwards towards the surface of the star
+ Calculates a few starting values at the stellar surface using an approximation where mass and luminosity are assumed constant and at their full values
+ Inputs these starting values into the same Adams-Bashforth integration algorithm and calculates inwards towards the centre of the star

There is some undeveloped code for calculations in transformed variables. This needs completion.

The code uses the following sources :

+ Figure 1 in Paper 3 gives the basic differential equations
+ Figure 2 in Paper 3 gives the constant density approximation equations for the stellar centre
+ Section 2.1.2 in Paper 3 gives the constant mass and luminosity approximation equations for the stellar surface
+ Section 2.1.3 in Paper 3 describes the Adams-Bashforth integration algorithm 

The code does not :

+ Attempt to join the outward and inward solutions to get a full solution
+ Do calculations for convective energy transfer
+ Do evolutionary calculations where the star develops in time and its chemical composition changes
+ Apply the Henyey method of stellar structure calculation
+ Do calculations with the stellar mass (instead of distance) as the independent variable.

The code was developed under VisualStudio v16.8.4. It is not guaranteed to be bug-free nor to be in the most elegant form.  
Please use and develop this code and the uploaded Papers 1-3 as you wish but credit me please : 

Peter Gardner : peter.gardner13@btopenworld.com, UK 07495 416625, also see my LinkedIn profile
