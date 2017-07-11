SUNDIALS: SUITE OF NONLINEAR AND DIFFERENTIAL/ALGEBRAIC EQUATION SOLVERS
------------------------------------------------------------------------


CVODE
=====

**CVODE** from LLNL is a solver for stiff and nonstiff ordinary differential equation
(ODE) systems (initial value problem) given in explicit form y’ = f(t,y).

This is a Github public reposoitory of LLNL's CVODE intended to make it easy to
add CVODE to existing projects as simply a git submodule. 

                            SUNDIALS 
    SUite of Nonlinear and DIfferential/ALgebraic equation Solvers
                   Release 2.6.2, August 2015
      Alan Hindmarsh, Daniel Reynolds, Radu Serban, Carol Woodward
           Center for Applied Scientific Computing, LLNL

The family of solvers referred to as SUNDIALS consists of the following solvers:
 ARKODE - for integration of ordinary differential equation systems (ODEs)
          ARKODE treats stiff, nonstiff and multi-rate ODE systems of the form
          y' = fe(t,y) + fi(t,y), y(t0) = y0
 CVODE  - for integration of ordinary differential equation systems (ODEs)
          CVODE treats stiff and nonstiff ODE systems of the form
          y' = f(t,y), y(t0) = y0
 CVODES - for integration and sensitivity analysis of ODEs
          CVODES treats stiff and nonstiff ODE systems of the form
          y' = f(t,y,p), y(t0) = y0(p)
 IDA    - for integration of differential-algebraic equation systems (DAEs)
          IDA treats DAE systems of the form
          F(t,y,y') = 0, y(t0) = y0, y'(t0) = y0'
 IDAS   - for integration and sensitivity analysis of DAEs
          IDAS treats DAE systems of the form
          F(t,y,y',p) = 0, y(t0) = y0(p), y'(t0) = y0'(p)
 KINSOL - for solution of nonlinear algebraic systems
          KINSOL treats nonlinear systems of the form
          F(u) = 0

The various solvers of this family share many subordinate modules.
For this reason, it is organized as a family, with a directory structure 
that exploits that sharing. Each individual solver includes documentation 
on installation, along with full usage documentation.

Warning to users who receive more than one of these individual solvers
at different times: The mixing of old and new versions of SUNDIALS may fail.  
To avoid such failures, obtain all desired solvers at the same time.

For installation directions see the file INSTALL_GUIDE.pdf.

Release history:


         
|  Date    | SUNDIALS release |  ARKODE  |   CVODE  | CVODES   |   IDA    |   IDAS   |  KINSOL  |
|----------|------------------|----------|----------|----------|----------|----------|----------|
| Sep 2016 |   2.7.0          |  1.1.0   |  2.9.0   |  2.9.0   |  2.9.0   |  1.3.0   |  2.9.0   |
| Aug 2015 |   2.6.2          |  1.0.2   |  2.8.2   |  2.8.2   |  2.8.2   |  1.2.2   |  2.8.2   |
| Mar 2015 |   2.6.1          |  1.0.1   |  2.8.1   |  2.8.1   |  2.8.1   |  1.2.1   |  2.8.1   |
| Mar 2015 |   2.6.0          |  1.0.0   |  2.8.0   |  2.8.0   |  2.8.0   |  1.2.0   |  2.8.0   |
| Mar 2012 |   2.5.0          |          |  2.7.0   |  2.7.0   |  2.7.0   |  1.1.0   |  2.7.0   |
| May 2009 |   2.4.0          |          |  2.6.0   |  2.6.0   |  2.6.0   |  1.0.0   |  2.6.0   |
| Nov 2006 |   2.3.0          |          |  2.5.0   |  2.5.0   |  2.5.0   |          |  2.5.0   |
| Mar 2006 |   2.2.0          |          |  2.4.0   |  2.4.0   |  2.4.0   |          |  2.4.0   |
| May 2005 |   2.1.1          |          |  2.3.0   |  2.3.0   |  2.3.0   |          |  2.3.0   |
| Apr 2005 |   2.1.0          |          |  2.3.0   |  2.2.0   |  2.3.0   |          |  2.3.0   |
| Mar 2005 |   2.0.2          |          |  2.2.2   |  2.1.2   |  2.2.2   |          |  2.2.2   |
| Jan 2005 |   2.0.1          |          |  2.2.1   |  2.1.1   |  2.2.1   |          |  2.2.1   |
| Dec 2004 |   2.0            |          |  2.2.0   |  2.1.0   |  2.2.0   |          |  2.2.0   |
| Jul 2002 |   1.0            |          |    2.0   |    1.0   |    2.0   |          |    2.0   |



The methods used in CVODE are variable-order, variable-step multistep
methods. For nonstiff problems, CVODE includes the Adams-Moulton formulas, with
the order varying between 1 and 12. For stiff problems, CVODE includes the
Backward Differentiation Formulas (BDFs) in so-called fixed-leading coefficient
form, with order varying between 1 and 5. For either choice of formula, the
resulting nonlinear system is solved (approximately) at each integration
step. For this, CVODE offers the choice of either functional iteration, suitable
only for nonstiff systems, and various versions of Newton iteration. In the
cases of a direct linear solver (dense or banded), the Newton iteration is a
Modified Newton iteration, in that the Jacobian is fixed (and usually out of
date). When using a Krylov method as the linear solver, the iteration is an
Inexact Newton iteration, using the current Jacobian (through matrix-free
products), in which the linear residual is nonzero but controlled.



The implicit nonlinear systems within implicit integrators are solved
approximately at each integration step using a modified Newton method, an
Inexact Newton method, or fixed-point solver (functional iteration). For the
Newton-based methods and the serial or threaded NVECTOR modules in SUNDIALS,
CVODE provides both direct (dense, band, or sparse) and preconditioned Krylov
iterative (GMRES, BiCGStab, TFQMR) linear solvers. When used with one of the
distributed parallel NVECTOR modules, including PETSc and hypre vectors, or a
user-provided vector data structure, only the Krylov solvers are available,
although a user may supply their own linear solver for any data structures if
desired.  For the serial vector structure, there is a banded preconditioner
module called CVBANDPRE for use with the Krylov solvers, while for the
distributed memory parallel structure there is a preconditioner module called
CVBBDPRE which provides a band-block-diagonal preconditioner.



For use with Fortran applications, a set of Fortran/C interface routines, called
FCVODE, is also supplied. These are written in C, but assume that the user
calling program and all user-supplied routines are in Fortran.


See the official Sundials
page <https://computation.llnl.gov/projects/sundials/cvode> 
