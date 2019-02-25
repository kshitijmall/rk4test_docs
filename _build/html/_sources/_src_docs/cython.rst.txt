Solution using Cython
=====================
1. EOM Function
----------------

.. code-block:: python

    cdef eom (double t, cnp.ndarray[cnp.float64_t,ndim=1] x):
    ################################################################################
    ## eom calculates the state derivative vectors at the given time and states
    #
    #  Author: Kshitij Mall
    #
    #  Modified: 24 February 2019
    #  Parameters:
    #
    #    Input: double t, the given time.
    #
    #           double vector x, the given state vector.
    #
    #    Output: double xdot: the state derivative vector
    #
    ################################################################################

      # Inputs to calculate the state derivatives
      cdef:
          double g = 9.80665 # Acceleration due to gravity, m/s^2
          double re = 6378000.0 # Radius of Earth, m
          double B = 157.0 # Ballistic Coefficient, kg/m^2
          double rho0 = 1.225 # Surface atmospheric density, kg/m^3
          double H = 7200.0 # Scale height, m
          cnp.npy_intp *dims = [3] # Dimension of the derivative vector
          cnp.ndarray[cnp.float64_t,ndim=1,mode='c'] xdot = \
            cnp.PyArray_EMPTY(1, dims, cnp.NPY_DOUBLE, 0)
      xdot[0] = x[1]*sin(x[2])
      xdot[1] = -rho0*exp(-x[0]/H)*pow(x[1],2)/(2*B)-g*sin(x[2])
      xdot[2] = (x[1]/(re+x[0]) - g/x[1])*cos(x[2])
      return xdot # Calculate the state derivative vector


2. RK4 File
------------

.. code-block:: python

    # cython: profile=True
    # cython: boundscheck=False
    # cython: wraparound=False
    # cython: cdivision=True
    # cython: profile=True
    cimport numpy as cnp

    from libc.math cimport (cos, sin, exp, pow) # Import all C-based libraries

    cnp.import_array() # This step is needed for ndarrays to work

    cpdef rk4 (double a, double b, int N, int M, cnp.ndarray[cnp.float64_t,ndim=1] alpha):
    ################################################################################
    ## RK4 calculates the state values for each time step
    #
    #  Author: Kshitij Mall
    #
    #  Modified: 24 February 2019
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

    # Declare and initialize all the necessary constants and variables for RK4 method
      cdef:
          double h = (b-a)/N # The step size
          cnp.npy_intp *dim1 = [N]
          cnp.npy_intp *dim2 = [M,N]
          cnp.ndarray[cnp.float64_t,ndim=1,mode='c'] t = \
           cnp.PyArray_EMPTY(1, dim1, cnp.NPY_DOUBLE, 0)
          cnp.ndarray[cnp.float64_t,ndim=2,mode='c'] w = \
           cnp.PyArray_EMPTY(2, dim2, cnp.NPY_DOUBLE, 0)
          Py_ssize_t i, j # Needed to handle 64-bit architectures
          cnp.ndarray[cnp.float64_t,ndim=1,mode='c'] k1, k2, k3, k4
          cnp.ndarray[cnp.float64_t,ndim=1,mode='c'] k2inp = \
           cnp.PyArray_EMPTY(1, dim1, cnp.NPY_DOUBLE, 0)
          cnp.ndarray[cnp.float64_t,ndim=1,mode='c'] k3inp = \
           cnp.PyArray_EMPTY(1, dim1, cnp.NPY_DOUBLE, 0)
          cnp.ndarray[cnp.float64_t,ndim=1,mode='c'] k4inp = \
           cnp.PyArray_EMPTY(1, dim1, cnp.NPY_DOUBLE, 0)

      t[0] = a # Initial time
      for j in range(M):
          w[j,0] = alpha[j] # Insert initial value input as the first row of the state
      # Obtain the states for the given times using rk4
      for i in range(N-1):
          k1 = eom(t[i], w[:,i]) # Calculate k1 for each time inputs
          for j in range(M):
              k2inp[j] = w[j,i]+0.5*k1[j]
          k2 = eom(t[i]+h/2, k2inp) # Calculate k2 for each time inputs
          for j in range(M):
              k3inp[j] = w[j,i]+0.5*k2[j]
          k3 = eom(t[i]+h/2, k3inp) # Calculate k3 for each time inputs
          for j in range(M):
              k4inp[j] = w[j,i]+k3[j]
          k4 = eom(t[i]+h, k4inp) # Calculate k4 for each time inputs
          t[i+1] = t[i] + h # Update the time values
          # Update the state matrix with k1, k2, k3, k4, and old state value
          for j in range(M):
              w[j,i+1] = w[j,i] + h*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0
      return t, w # Return the time and state vectors

3. Main Input File
-------------------

.. code-block:: python

    import sys
    import time
    import numpy as np
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
        print ( '  Test the RK4 Function with Cython.' )
        print ( '' )
        # Write the necessary inputs
        vatm = 11060 # Entry Velocity, m/s
        hatm = 80000 # Entry Height, m
        gamma0 = -50/180*np.pi # Initial flight path angle, rad
        t0 = 0 # Initial time
        tf = 212.2 # Final time
        step = 1000 # Time steps
        S = 3 # Number of states
        init = np.array([hatm,vatm,gamma0],np.float64) # Initial state vector

        try:
            tic = time.time() # Start the timer
            # Obtain the states for the given times using rk4
            [tout,y] = rk4.rk4(t0, tf, step, S, init)
            elapsed = time.time() - tic # Calculate the elapsed time
            # Print the computation time
            print('Time taken by python with cython code:',elapsed)
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

4. Setup File
--------------

.. code-block:: python

    # distutils: define_macros=CYTHON_TRACE_NOGIL=1
    from distutils.core import setup
    from distutils.extension import Extension
    from Cython.Build import cythonize
    import numpy as np
    import Cython

    #This will generate HTML to show where there are still pythonic bits hiding out
    Cython.Compiler.Options.annotate = True

    ext_modules = [
        Extension("rk4new",
                  sources=["rk4.pyx"],
                  libraries=["m"]  # Unix-like specific
                  )
    ]

    setup(name="rk4project",
          ext_modules=cythonize(['rk4.pyx']),
          include_dirs=[np.get_include()]
    )
