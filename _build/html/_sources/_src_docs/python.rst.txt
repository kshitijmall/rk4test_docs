Solution using Python
=====================
1. EOM File
------------

.. code-block:: python

    # Import the necessary packages
    import numpy as np
    def eom (t,x):
    ################################################################################
    ## eom calculates the state derivative vectors at the given time and states
    #
    #  Author: Kshitij Mall
    #
    #  Modified: 18 February 2019
    #  Parameters:
    #
    #    Input: double t, the given time.
    #
    #           double x, the given state vector.
    #
    #    Output: double xdot: the state derivative vector
    #
    ################################################################################
        # Inputs to calculate the state derivatives
        g = 9.80665
        re = 6378000
        B = 157
        rho0 = 1.225
        H = 7200
        xdot = [x[1]*np.sin(x[2]),\
               -rho0*np.exp(-x[0]/H)*x[1]**2/(2*B)-g*np.sin(x[2]),\
               (x[1]/(re+x[0]) - g/x[1])*np.cos(x[2])]
        return xdot # Calculate the state derivative vector

2. RK4 File
------------

.. code-block:: python

    # Import the necessary packages
    import numpy as np
    import eom
    def rk4 (a, b, N, M, alpha):
    ################################################################################
    ## RK4 calculates the state values for each time step
    #
    #  Author: Kshitij Mall
    #
    #  Modified: 18 February 2019
    #  Parameters:
    #
    #    Input: double a, the initial time.
    #
    #           double b, the final time.
    #
    #           integer N, the time step.
    #
    #           integer M, the number of states.
    #
    #           double alpha, the initial state vector.
    #
    #    Output: double t, w: the fourth-order Runge-Kutta solution
    #
    ################################################################################
        h = float((b-a))/N # The step size
        t = np.empty ( N ) # Initialize the time vector with zero values
        t[0] = a # Initial time
        w = np.empty ( (M, N) ) # Initialize the state matrix with zero values
        w[:,0] = alpha # Insert initial value input as the first row of the state
        # Obtain the states for the given times using rk4
        for i in range(N-1): # Calculate k1, k2, k3, and k4 for each time inputs
            k1 = h*np.ones (M)*eom.eom(t[i], w[:,i])
            k2 = h*np.ones (M)*eom.eom(t[i]+h/2, w[:,i]+0.5*k1)
            k3 = h*np.ones (M)*eom.eom(t[i]+h/2, w[:,i]+0.5*k2)
            k4 = h*np.ones (M)*eom.eom(t[i]+h, w[:,i]+k3)\
            # Update the state matrix with k1, k2, k3, k4, and old state value
            w[:,i+1] = w[:,i] + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
            t[i+1] = t[i] + h # Update the time values
        return t, w # Return the time and state vectors

3. Main Input File
-------------------

.. code-block:: python

    # Import the necessary packages
    import sys
    import time
    from math import pi
    import platform
    import rk4
    def input():
    ################################################################################
    # main
    #
    #  Author: Kshitij Mall
    #
    #  Modified: 24 February 2019
    #  Parameters:
    #
    #    Input: None
    #
    #    Output: real [t,y], the fourth-order Runge-Kutta solution
    #
    ################################################################################
        print ( '' )
        print ( '  Python version: %s' % ( platform.python_version ( ) ) )
        print ( '  Test the RK4 Function.' )
        print ( '' )
        # Write the necessary inputs
        vatm = 11060 # Entry Velocity, m/s
        hatm = 80000 # Entry Height, m
        gamma0 = -50/180*pi # Initial flight path angle, rad
        t0 = 0 # Initial time
        tf = 212.2 # Final time
        step = 1000 # Time steps
        S = 3 # Number of states
        init = [hatm,vatm,gamma0] # Initial state vector

        try:
            tic = time.time() # Start the timer
            # Obtain the states for the given times using rk4
            [t,y] = rk4.rk4(t0, tf, step, S, init)
            elapsed = time.time() - tic # Calculate the elapsed time
            # Print the computation time
            print('Time taken by pure python code:',elapsed)
        except:
            # In case of an unexpected error catch and raise an exception
            print("Unexpected error:", sys.exc_info()[0])
            raise

    """The following if condition allows this python module to be imported by other modules
    without calling main function. If desired, this main function can be called by
    the other module that imports this module.
    """
    if ( __name__ == '__main__' ):
       input ( )
