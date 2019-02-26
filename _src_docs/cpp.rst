Solution using C++
==================
1. EOM Function
---------------

.. code-block:: c

    //******************************************************************************

    double * eom(double t, double x[]) {

    //******************************************************************************
    //
    //  Author:
    //
    //    Kshitij Mall
    //
    //  Modified:
    //
    //    18 February 2019
    //
    //  Purpose:
    //
    //    eom is the function file to obtain the derivative vector at a given time.
    //
    //  Inputs:
    //
    //    t: Time at which derivative needs to be evaluated
    //    x[]: starting state values for h, v, and fpa at the input time
    //
    //  Output:
    //
    //    xdot: Vector output with derivatives for each of the three input states
    //
    //******************************************************************************
    	const double g = 9.80665; // Acceleration due to gravity, m/s^2
    	const double re = 6378000.0; // Radius of Earth, m
    	const double B = 157.0; // Ballistic coefficient, kg/m^2
    	const double rho0 = 1.225; // Surface atmospheric density, kg/m^3
    	const double H = 7200.0; // Scale height, m

      static double xdot[S]; // Declare the derivative vector

      xdot[0] = x[1]*sin(x[2]); // hdot
      xdot[1] = -rho0*exp(-x[0]/H)*pow(x[1],2.0)/(2.0*B) - g*sin(x[2]); // vdot
      xdot[2] = (x[1]/(re+x[0]) - g/x[1])*cos(x[2]); // fpa dot

      return xdot;
    }

2. RK4 Function
----------------

.. code-block:: c

    //******************************************************************************

    double ** rk4(double a, double b, double alpha[]){

    //******************************************************************************
    //
    //  Author:
    //
    //    Kshitij Mall
    //
    //  Modified:
    //
    //    18 February 2019
    //
    //  Purpose:
    //
    //    rk4 is the function file to obtain the state vector at all input times.
    //    A 4th order Runge Kutta method is used in this file.
    //
    //  Inputs:
    //
    //    a: Initial time
    //		b: Final time
    //    alpha[]: Initial state vector
    //
    //  Output:
    //
    //    state: A matrix of the state outputs at different times
    //
    //******************************************************************************
      double ** state = new double*[N]; // Create N pointers to the rows of state

      double h = double((b-a)/N); // The step size
    	// Initialize useful variables for the Runge Kutta formula
      double * k1, * k2, * k3, * k4, t[N+1], k2inp[S], k3inp[S], k4inp[S];

      for (int i = 0; i < N+1; ++i) {
           state[i] = new double[N]; // Insert columns of states for each time input
           if(i == 0){
    				 // Initialize the time and state arrays
             t[i] = a;
    				 for (int j = 0; j < S; ++j) {
             state[i][j] = alpha[j];
    			 	}
           }
           else{
             // Determine k1
             k1 = eom(t[i-1], state[i-1]);

             // Determine k2
             for (int j = 0; j < S; ++j) {
              k2inp[j]  = state[i-1][j] + 0.5*k1[j];
             }
             k2 = eom(t[i-1] + h/2.0, k2inp);

             // Determine k3
             for (int j = 0; j < S; ++j) {
              k3inp[j]  = state[i-1][j] + 0.5*k2[j];
             }
             k3 = eom(t[i-1] + h/2.0, k3inp);

             // Determine k4
             for (int j = 0; j < S; ++j) {
              k4inp[j]  = state[i-1][j] + k3[j];
             }
             k4 = eom(t[i-1] + h, k4inp);

             // Determine the new state from k1, k2, k3, k4, and the old state
             for (int j = 0; j < S; ++j) {
               state[i][j] = state[i-1][j] + h*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0;
             }

    				 t[i] = t[i-1] + h; // Determine the new time
           }
      }

      return state; // Return the state vector
    }


3. printResult Function
-----------------------

.. code-block:: c

    //******************************************************************************

    void printResult(double ** state){

    //******************************************************************************
    //
    //  Author:
    //
    //    Kshitij Mall
    //
    //  Modified:
    //
    //    18 February 2019
    //
    //  Purpose:
    //
    //    printResult is the function file to obtain the derivative vector at a given time.
    //
    //  Inputs:
    //
    //    state: State result vector obtained from rk4 function
    //
    //  Output:
    //
    //    Returns void
    //
    //******************************************************************************
      for (int i = 0; i < N+1; ++i) {
        for (int j = 0; j < S; ++j) {
    			// Print all the state vector values at all input times
          cout << "state[" <<i<<"][" <<j <<"] = " << state[i][j] << " ";
        }
        cout << endl;
      }
    }


4. Main File
-------------

.. code-block:: c

    #include <iostream>
    // To use mathematical quantities like PI we need the following two lines
    #define _USE_MATH_DEFINES
    #include <math.h>
    #include <time.h> // This header file is needed to calculate computation time

    using namespace std; // This helps to avoid writing std:: at several places

    // Define the constants used for the hypersonic problem
    const int N = 1000; // Time steps
    const int S = 3; // Number of states

    // Declare the functions used inside main
    double * eom(double t, double x[]);
    double ** rk4(double a, double b, double alpha[]);
    void printResult(double ** result);

    //******************************************************************************

    int main (int argc, const char * argv[]) {

    //******************************************************************************
    //
    //  Author:
    //
    //    Kshitij Mall
    //
    //  Modified:
    //
    //    24 February 2019
    //
    //  Purpose:
    //
    //    main is the main program for this project.
    //
    //  Inputs:
    //
    //    argc: Number of command line arguments
    //    argv: value of array of command line arguments
    //
    //  Output:
    //
    //    Returns 0
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license
    //
    //******************************************************************************
    cout << "\n" << "  C++ version\n";

    const double vatm = 11060.0; // Entry Velocity, m/s
    const double hatm = 80000.0; // Entry Height, m
    const double gamma0 = -50.0/180.0*M_PI; // Initial flight path angle, rad

    const double t0 = 0.0; // Initial time
    const double tf = 212.2; // Final time

    double init[] = {hatm,vatm,gamma0}; // Initial guess

    // Get starting timepoint
    clock_t t;
    t = clock();

    // Call the RK4 function to obtain the ODE solution
    double ** result = rk4(t0, tf, init);

    // Get ending timepoint
    t = clock() - t;

    // Display the computation time
    cout << "  Computation Time: " << t*1.0/CLOCKS_PER_SEC << " seconds" << endl;

    //printResult(result); // Print the state vector
    return 0;
    }
