.. _schroedinger-poisson:

Schrödinger-Poisson Solver
==========================

.. contents::
   :depth: 2
   :local:

Simulation Constants
--------------------

.. automodule:: solver.schroedinger_poisson.common

.. autoclass:: sim_constants
   :show-inheritance:
   :members:

Solvers
-------

.. automodule:: solver.schroedinger_poisson.solver

Transfer-Matrix Solver
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: tm_solver
   :show-inheritance:
   :members:

Poisson solver
^^^^^^^^^^^^^^

.. autoclass:: poisson_solver
   :show-inheritance:
   :members:


Boundary Conditions
-------------------

.. automodule:: solver.schroedinger_poisson.common

.. autoclass:: boundary_condition
   :show-inheritance:
   :members:

.. autoclass:: ez_base_transform
   :show-inheritance:
   :members:

.. autoclass:: tb_base_transform
   :show-inheritance:
   :members:

.. automodule:: solver.schroedinger_poisson.boundary_conditions

.. autoclass:: boundary_condition_ext
   :show-inheritance:
   :members:

.. autoclass:: boundary_condition_ez
   :show-inheritance:
   :members:

.. autoclass:: boundary_condition_QCD
   :show-inheritance:
   :members:

.. autoclass:: boundary_condition_tb
   :show-inheritance:
   :members:
