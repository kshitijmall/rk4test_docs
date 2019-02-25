.. RK4Test documentation master file, created by
   sphinx-quickstart on Mon Feb 25 03:36:50 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

RK4 Solution for a Hypersonics Initial Value Problem
==========================================================
In this mini-project an initial value problem (IVP) from hypersonics domain is
integrated using the RK4 method. The equations of motion (EOMs) for a ballistic
hypersonic vehicle reentering the Earth used in this IVP are as follows.

.. math ::
  \dot{h} &= v\sin{\gamma}

  \dot{v} &= \dfrac{-\rho_{0}e^{\frac{-h}{H}}v^{2}}{2B}-g\sin{\gamma}

  \dot{\gamma} &= \left(\dfrac{v}{re + h} - \dfrac{g}{v}\right)\cos{\gamma}

where
:math:`h` is the altitude,
:math:`v` is the velocity,
:math:`\gamma` is the flight path angle,
:math:`\rho_{0}` is the surface atmospheric density,
:math:`H` is the scale height,
:math:`re` is the radius of Earth,
:math:`B` is the ballistic coefficient of the hypersonic vehicle
and
:math:`g` is the acceleration due to Earth's gravity.

The initial conditions for this IVP are as follows.

.. math ::
  h(t_{0}) &= 80000\ m

  v(t_{0}) &= 11060\ \dfrac{m}{s}

  \gamma(t_{0}) &= -50\ deg

The values of the constants used in this IVP are as follows.

.. math ::
  g &= 9.80665\ \dfrac{m}{s^2}

  re &= 6378000\ m

  B &= 157\ \dfrac{kg}{m^2}

  \rho_{0} &= 1.225\ \dfrac{kg}{m^3}

  H &= 7200\ m

The EOMs are integrated from an initial time of 0 s to a terminal time of 212.2
s using the following RK4 method with 1000 steps.

.. math ::
  k_{1} &= h~\text{eom}(t(i),w(:,i))

  k_{2} &= h~\text{eom}\left(t(i) + \frac{h}{2}, w(:,i) + 0.5~k_{1}\right)

  k_{3} &= h~\text{eom}\left(t(i) + \frac{h}{2}, w(:,i) + 0.5~k_{2}\right)

  k_{4} &= h~\text{eom}(t(i) + h, w(:,i) + k_{3})

  w(:,i+1) &= w(:,i) + \frac{k_{1} + 2k_{2} + 2k_{3} + k_{4}}{6}

  t(i+1) &= t(i) + h

where
:math:`h` is the step size,
:math:`t(i)` is the time at the nth step,
:math:`w(:,i)` is the state vector at the nth step,
eom is the function that calculates the state vector derivatives,
:math:`k_{1}`, :math:`k_{2}`, :math:`k_{3}`, and :math:`k_{4}` are the vectors used in the RK4 method.

The resulting energy plot for this IVP using RK4 method is shown below.

.. image:: ./energy_result.png
  :width: 800
  :alt: Alternative text

Documentation
-------------

.. toctree::
   :maxdepth: 1
   :titlesonly:

   ./_src_docs/python.rst
   ./_src_docs/cython.rst
   ./_src_docs/cpp.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
