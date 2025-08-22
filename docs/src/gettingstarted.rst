.. _gettingstarted:

Getting started
===============

This tutorial outlines the first steps with monaco. The first task is usually
the setup of a new simulation such as the transport simulation of a certain new
structure.

Overview of the simulation tools
--------------------------------

The Schrödinger-Poisson solver and the self-consistent Rate-equation solver are
the two main components of such a transport simulation. Both programs are
written in MATLAB in an object-oriented programming style.

Implementations
^^^^^^^^^^^^^^^

- **Schrödinger-Poisson:**
  The Schrödinger equation is solved using the transfer-matrix method, with
  either extended or tight-binding basis as boundary conditions for the
  wavefunctions.
  To include the effect of space charges and free carriers, the Schrödinger
  equation is coupled to the Poisson equation in a self-consistent manner.

  Path: ``solver/schroedinger_poisson``

- **Rate-equations:**
  The transport through the device is simulated based on a rate-equation
  description, where the confined energy states act as laser levels.
  The scattering rates between these levels are calculated from the wave
  functions and material parameters using Fermi's golden rule.
  The following scattering mechanisms are implemented and can be
  enabled/disabled in a simulation:
  LO-Phonon, acoustic phonon, alloy, impurity, interface roughness and
  electron-electron scattering, as well as, tunneling, stimulated emission and
  absorption.

  Path: ``solver/rate_equation``

- **Self-consistent scheme**:
  The solver class ``SP_Rate_solver`` offers a convenient way to simulate the
  charge carrier transport through the device in a self-consistent manner.
  Here, the Schrödinger-Poisson equation and the rate equations are solved
  iteratively until convergence is reached.

 Path: ``solver/transport``

Restrictions
^^^^^^^^^^^^

The project is in active development. While we aim to make the code as general
and flexible as possible, the following restriction apply:

.. note::

  The MATLAB code has not been tested with Octave.
  In particular, the object-oriented part is not likely to work with Octave.

How to setup a simulation
-------------------------

First, find a place to store your project files. This could be in
`setups/my_little_project` or any place outside the source directory (in the
latter case make sure that MATLAB will find the framework in `setups/matlab`).

Then, we recommend that you copy an existing yaml-file defining a specific setup
and make the necessary changes. Describe the design of the active region by
setting up the materials and adding layers to the structure `period`.
Using this information, you can create a `device`.

Every `device` can be simulated with a certain `scenario`, which defines the
environment such as the operating temperature and the applied bias.

Finally, solver-specific settings can be specified.

.. note::

  - Make sure that the period starts with the highest barrier.

  - Be cautious when using more than two different materials as this has not
    been tested thoroughly.

  - Use at least 5 periods for the device structure.

  - Setting ``num_wavefct`` specifies the number of wave functions that are
    considered by the solvers. You have to inspect the results of the
    Schrödinger-Poisson solver first to determine how many wave functions seem
    relevant.

How to run the simulation
-------------------------------

First, read the yaml input file, which defines the ``device`` and ``scenario``.

.. code::

  [d, s] = read_input_file("benz2007");

Then, run the Schrödinger-Poisson solver. At this point, you can inspect the
wave functions and verify the device description in your yaml file.

.. code::

  SP_solver = tm_solver(d, s);
  [eigen, cond_band] = SP_solver.solve();
  eigen.plot_wavefunctions(cond_band);

Now you can initialize the carrier transport solver by choosing which
scattering mechanisms should be considered in the simulation. For simulations
that take the electromagnetic field in the cavity into account, include
``optical field`` in the list of scattering mechanisms. After initialization,
you can start the the self-consistent loop.

.. code::

  scattering_mech = {"tunneling", "alloy disorder", "impurity", ...
    "acoustic phonon", "lo phonon emission", "lo phonon absorption", ...
    "interface roughness", "electron electron"};
  SPR_Solver = SP_Rate_solver(d, s, scattering_mech);

  [cond, eigen, J_new, dist, deph, sc, gain, rateSolver, m] = SPR_Solver.solve();

After the simulation has finished, the results are returned as postprocessing
classes. Use the ``hdf5_write()`` function to save the simulation results in
hdf5 file format.

.. code::

  filename = 'simulation_results.h5'
  hdf5_write(filename, cond, eigen, dist, deph, sc)
