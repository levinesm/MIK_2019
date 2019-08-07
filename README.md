# MIK_2019
Modified Inverse Kirkendall (MIK) Code

:::General Program Information::::

Program Title: Main (local.f95)

Programming Language: FORTRAN 95

Subroutines:
  Initialize – provides formats for all inputs. Reads in values from input file. Begins writing output file (writes initial values to file). Calculates binding energies and migration energies. Sets up mesh? Calculates diffusion coefficients, thermal concentration of vacancies.
  Prep – sets up DLSODE solver
  Output – solves for concentration of A atoms, B atoms, C atoms, vacancies, and interstitials. Writes outputs to Unit 6.
  FEX – computes derivatives using DLSODE solver
  JEX – required as a dummy variable by the DLSODE solver. Given MITER=2, the analytical Jacobian matrix is generated internally. Subroutine is required for the solver to run, but it does not have any functionality.

Other files:
  Unit 5 – perks.in
  Unit 6 – perks.out
  Unit 8 – perks.err

Variables:
  Mesh Definition:
    R1 – upper bound of region 1 (gb to distance 1 (nm))
    R2 – upper bound of region 2 (distance 1 to distance 2 (nm))
    R3 – upper bound of region 3 (distance 2 to distance 3 (nm))
    N1 – number of steps (subdivisions) for mesh in region 1
    N2 – number of steps (subdivisions) for mesh in region 2
    N3 – number of steps (subdivisions) for mesh in region 3
  Tolerances:
    EPSR – local relative error tolerance for the solution (in DLSODE, RTOL)
    EPSA – local absolute error tolerance for the solution (in DLSODE, ATOL)
  Irradiation:
    DISPRT – displacement rate (dpa/s)
    ETAV – displacement efficiency for production for freely migrating vacancies
    ETAI – displacement efficiency for production for freely migrating interstitials
    DOSE – dose (dpa)
    TEMPC – irradiation temperature (°C)
  Concentrations:
    CONCB – concentration of B atoms
    CONCC – concentration of C atoms
  Diffusion:
    NAT – atom number density
    LAMBDA – jump distance
    FAV – jump correlation factor of A atoms to vacancies
    FBV – jump correlation factor of B atoms to vacancies
    FCV – jump correlation factor of C atoms to vacancies
    FI – atom-interstation jump correlation factor
    WAV – relative vacancy jump frequency for A atoms (normalized to C atoms)
    WBV – relative vacancy jump frequency for B atoms (normalized to C atoms)
    WCV – relative vacancy jump frequency for C atoms (normalized to C atoms)
    WAI – relative interstitial jump frequency for A atoms (normalized to C atoms)
    WBI – relative interstitial jump frequency for B atoms (normalized to C atoms)
    WCI – relative interstitial jump frequency for C atoms (normalized to C atoms)
    NUOV – Debye frequency for vacancies
    NUOI – Debye frequency for interstitials
  Cohesive energies:
    ECOHA – cohesive energy of A atoms (when A=Fe, includes free energy difference between FCC and BCC)
    ECOHB – cohesive energy of B atoms (when B=Cr, includes free energy difference between FCC and BCC)
    ECOHC – cohesive energy of C atoms
  Migration energies:
    EMIA – migration energy of A atoms by interstitials
    EMIB – migration energy of B atoms by interstitials
    EMIC – migration energy of C atoms by interstitials
    SV – entropy for vacancy formation
    EMA – migration energy of A atoms by vacancies
    EMB – migration energy of B atoms by vacancies
    EMC – migration energy of C atoms by vacancies
  Formation energies:
    EFA – vacancy formation energy for A atoms
    EFB – vacancy formation energy for B atoms
    EFC – vacancy formation energy for C atoms
    EFGB – vacancy formation energy for grain boundary atoms
  Ordering energies:
    EORDAB – ordering energy for A-B atom pairs
    EORDAC – ordering energy for A-C atom pairs
    EORDBC – ordering energy for B-C atom pairs
  Material Properties:
    BIASV – bias factor for vacancies
    BIASI – bias factor for interstitials
    Z – recombination volume
    AL – thermodynamic factor
    DISL – dislocation density
  Program variables:
    PSTOP – program stop
    NOUT – number of output times
    TOUTPT – time outputs
    ISTEP – counts iterations
    IP – an arbitrary number that determines the size of all the arrays
