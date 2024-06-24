// Copyright (C) 2021 Yin Shi
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

// No need to include mpi.h, because dolfin.h already includes it
#include <math.h>
#include <fstream>
#include <dolfin.h>
#include "imt.h"
#include "pugixml.hpp"
#include "input.h"
#include <mpi.h>

using namespace dolfin;

// Use operator double() to cast dolfin Constant object to doulbe number

// Define a convenient function to split a multicomponent Function either by shallow copy
// or deep copy. Returns a shared pointer of i th subfunction. If deep copying, the copying
// memory is allocated by this function and its shared_prt is returned. This function can 
// be directly copied to other dolfin projects to use.
std::shared_ptr<Function> split(const Function& u, std::size_t i, bool deepcopy=false)
{
  auto usub_shallow = std::make_shared<Function>(u, i);  // Shallow copy, whose content cannot be modified inherently
  if (deepcopy)
  {
    // Deep copy, standalone, whose content can be modified. This shallow-then-deep copying can 
    // actually be used to create Functions on sub function space.
    auto usub_deep = std::make_shared<Function>(*usub_shallow);  
    return usub_deep;
  }
  return usub_shallow;
}

// Initial conditions
class InitialConditions : public Expression
{
  private:
  double _eta_i, _mu_i, _gamma_ei, _gamma_hi, _phi_i, _Ly, _T_i, _curr_i;

  public:
  InitialConditions(double eta_i, double mu_i, double gamma_ei, double gamma_hi, 
                    double phi_i, double Ly, double T_i, double curr_i) : Expression(7)
  {
    _eta_i = eta_i;
    _mu_i = mu_i;
    _gamma_ei = gamma_ei;
    _gamma_hi = gamma_hi;
    _phi_i = phi_i;
    _Ly = Ly;
    _T_i = T_i;
    _curr_i = curr_i;
  }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = _eta_i;
    values[1] = _mu_i;
    values[2] = _gamma_ei;
    values[3] = _gamma_hi;
    values[4] = _phi_i - _phi_i*x[1]/_Ly;
    values[5] = _T_i;
    values[6] = _curr_i;
  }
};

// Define Tc variation field
class TcVariation : public Expression
{
  private:
  double _TCIMP0, _RIMPDIS, _DWWIDTH, _Lx;
  
  public:
  TcVariation(double TCIMP0, double RIMPDIS, double DWWIDTH, double Lx) : Expression()  // Initialize single valued Expression
  {
    _TCIMP0 = TCIMP0;
    _RIMPDIS = RIMPDIS;
    _DWWIDTH = DWWIDTH;
    _Lx = Lx;
  }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = _TCIMP0*(-tanh(2*(sqrt(pow(x[0] - _Lx/2, 2) + pow(x[1], 2)) - _RIMPDIS)/_DWWIDTH) + 1)/2;
  }
};

// Define boundaries
class BoundaryY : public SubDomain
{
  public:
  BoundaryY(double y) : SubDomain()
  {
    _y = y;
  }

  private:
  double _y;

  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary && near(x[1], _y);
  }
};


// Define function for ramping up applied voltage
double ramped_delV(double t, double tramp, double delV_i, double delV0)
{
  double _delV;
  if (t <= tramp)
  {
    _delV = delV_i + (delV0 - delV_i) * t / tramp;
  }
  else
  {
    _delV = delV0;
  }
  return _delV;
}

// Calculate relative error between two solutions on the same function space using specified norm type 
double rel_err(std::shared_ptr<const Function> u1, std::shared_ptr<const Function> u2, imt::Form_scalernorm& scalernorm)
{
  Function delu(*u1);  // Copy constructor
  std::shared_ptr<const Function> s(nullptr);
  double a, b, r=0.0;

  *(delu.vector()) -= *(u2->vector());  // This is the difference Function between u1 and u2

  for (std::size_t i=0; i<6; ++i)
  {
    s = split(delu, i);
    scalernorm.s = s;
    a = assemble(scalernorm);
    s = split(*u1, i);
    scalernorm.s = s;
    b = assemble(scalernorm);
    r += (a / (b + 1e-15));
  }

  r = sqrt(r / 6);
  
  return r;
}

void save_sol(std::vector<File*> files, std::shared_ptr<const Function> u, const Function& ca, const double units[6], double t)
{
  std::shared_ptr<Function> tmp(nullptr);

  tmp = split(*u, 0);
  tmp->rename("eta", "eta");
  *(files[0]) << std::pair<const Function*, double>(&(*tmp), t);

  tmp = split(*u, 1);
  tmp->rename("mu", "mu");
  *(files[1]) << std::pair<const Function*, double>(&(*tmp), t);

  /*
  tmp = split(*u, 2);
  tmp->rename("gamma_e", "gamma_e");
  *(files[2]) << std::pair<const Function*, double>(&(*tmp), t);
  tmp = split(*u, 3);
  tmp->rename("gamma_h", "gamma_h");
  *(files[3]) << std::pair<const Function*, double>(&(*tmp), t);
  */

  tmp = split(*u, 4, true);  // Deep copy
  *(tmp->vector()) *= units[4];
  tmp->rename("phi", "phi");
  *(files[2]) << std::pair<const Function*, double>(&(*tmp), t);

  tmp = split(*u, 5, true);  // Deep copy
  *(tmp->vector()) *= units[5];
  tmp->rename("T", "T");
  *(files[3]) << std::pair<const Function*, double>(&(*tmp), t);
  
  *(files[4]) << std::pair<const Function*, double>(&(ca), t);
}



int main()
{
  MPI_Init(NULL, NULL);

  std::clock_t c_start = std::clock();  // Record calculation start time
  int worldsize, rank;

  MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "User message ===> Processor reporting: " << rank << "/" << worldsize << std::endl;

  // Define file related variables
  pugi::xml_document doc;
  pugi::xml_node docroot;
  std::string figdir = "solution/";

  //-----------------------------------------------------------------------
  // Define units (in SI units) to renormalize quantities
  //-----------------------------------------------------------------------
  const double LUNIT = 1e-9;
  const double TUNIT = 1e-9;
  const double TEMPUNIT = 338.0;
  const double EUNIT = 1.3806504e-23 * TEMPUNIT;
  const double VUNIT = 1e-3;
  const double CUNIT = EUNIT/VUNIT;
  const double MUNIT = EUNIT*pow(TUNIT/LUNIT, 2);
  const double RUNIT = VUNIT/(CUNIT/TUNIT);
  const double UNITS[6] = {1.0, 1.0, 1.0, 1.0, VUNIT, TEMPUNIT}; // Collect units for the 6 solution Functions

  if (!doc.load_file("input.xml")) 
  {
    std::cout << "Error: Cannot load file input.xml" << std::endl;
    return -1;
  }
  docroot = doc.child("input");
  //----------------------------------------------------------------
  // Define computational parameter Constant
  //----------------------------------------------------------------
  auto dt = std::make_shared<Constant>(1e-6);
  auto Lx = std::make_shared<Constant>(50e-9/LUNIT);
  readxml_bcast(*Lx, docroot, "external/Lx", MPI_COMM_WORLD, rank, LUNIT/1e-9);
  auto Lz = std::make_shared<Constant>(36e-9/LUNIT);
  readxml_bcast(*Lz, docroot, "external/Lz", MPI_COMM_WORLD, rank, LUNIT/1e-9);
  
  //----------------------------------------------------------------
  // Define physical property parameter Constant
  //----------------------------------------------------------------
  auto KB = std::make_shared<Constant>(1.3806504e-23/EUNIT*TEMPUNIT);
  auto ECHARGE = std::make_shared<Constant>(1.602176487e-19/CUNIT);
  const double UCVOL = 59e-30/pow(LUNIT, 3);
  auto TC = std::make_shared<Constant>(338.0/TEMPUNIT);
  auto T1 =  std::make_shared<Constant>(275.0/TEMPUNIT);   // 0.81361
  auto AN1 = std::make_shared<Constant>(2.05714*double(*KB)*double(*TC)/UCVOL);  // 34.867
  auto AN2 = std::make_shared<Constant>((-0.623108 + 0.121228)/2*double(*KB)*double(*TC)/UCVOL);
  auto AN3 = std::make_shared<Constant>((0.330568 + 4.18947)/4*double(*KB)*double(*TC)/UCVOL);
  auto T2 = std::make_shared<Constant>(270.0/TEMPUNIT);
  auto AU1 = std::make_shared<Constant>(3.94286*double(*KB)*double(*TC)/UCVOL);
  auto AU2 = std::make_shared<Constant>((1.36767 - 3.67915)/2*double(*KB)*double(*TC)/UCVOL);
  auto AU3 = std::make_shared<Constant>((0.4 + 2)/4*double(*KB)*double(*TC)/UCVOL);
  auto GNU1 = std::make_shared<Constant>(0.3*double(*KB)*double(*TC)/UCVOL);
  auto GNU2 = std::make_shared<Constant>((0.2 - 1.5 + 0.3/2)/2*double(*KB)*double(*TC)/UCVOL);
  auto GNU3 = std::make_shared<Constant>((0.05 + 2.0)/2*double(*KB)*double(*TC)/UCVOL);
  auto CHI = std::make_shared<Constant>(0.286*1.602176487e-19/EUNIT);
  readxml_bcast(*CHI, docroot, "internal/gapcoeff", MPI_COMM_WORLD, rank, EUNIT/1.602176487e-19);
  auto KN = std::make_shared<Constant>(1.0/((1e-12/TUNIT)*double(*AN1)*(double(*TC)-double(*T1))/double(*TC)));
  auto KAPPAN = std::make_shared<Constant>(1.0*(1.602176487e-19/1e-9)/EUNIT*LUNIT);
  auto KU = std::make_shared<Constant>(1.0/((10e-15/TUNIT)*double(*CHI)*2*(0.16/UCVOL)));
  auto KAPPAU = std::make_shared<Constant>(1.0*(1.602176487e-19/1e-9)/EUNIT*LUNIT);
  auto NC = std::make_shared<Constant>(2*pow( 65*(9.10938215e-31/MUNIT) * double(*KB)*double(*TC) / (2*M_PI * pow(1.054571628e-34/EUNIT/TUNIT, 2)), 1.5 ));
  auto NV = std::make_shared<Constant>(double(*NC));
  auto MEC = std::make_shared<Constant>(5e-5/pow(LUNIT, 2)*VUNIT*TUNIT);
  auto MEA = std::make_shared<Constant>(double(*MEC)*0.5);
  auto MHC = std::make_shared<Constant>(double(*MEC)/1.2);
  auto MHA = std::make_shared<Constant>(double(*MHC)*0.5);
  auto KEH0 = std::make_shared<Constant>(1.0/( 2 * sqrt(double(*NC)*double(*NV))*exp(-double(*CHI)*pow(0.827987, 2)/(double(*KB)*322.0/TEMPUNIT)) * (14.235e-6/TUNIT) ));
  readxml_bcast(*KEH0, docroot, "internal/ehrecrate", MPI_COMM_WORLD, rank, pow(LUNIT, 3)/TUNIT/(UCVOL/1e-9));
  auto PERMITTIVITY = std::make_shared<Constant>(60*(8.854187817e-12*VUNIT/CUNIT*LUNIT));
  auto CPV = std::make_shared<Constant>(690*4340/EUNIT*pow(LUNIT, 3)*TEMPUNIT);
  auto THETA = std::make_shared<Constant>(6/EUNIT*TUNIT*LUNIT*TEMPUNIT);
  auto HTRAN = std::make_shared<Constant>(3e6/EUNIT*TUNIT*pow(LUNIT, 2)*TEMPUNIT);
  readxml_bcast(*HTRAN, docroot, "external/heatdiss", MPI_COMM_WORLD, rank, EUNIT/(TUNIT*pow(LUNIT, 2)*TEMPUNIT));
  auto CHP_IN = std::make_shared<Constant>(0.0);   // Intrinsic chemical potential
  
  //----------------------------------------------------------------
  // Define computational parameters
  //----------------------------------------------------------------
  int nx = 50, ny = 20, max_div = 100, max_iter = 30, rseed = 11793, kry_max_iter = 1000;
  readxml_bcast(nx, docroot, "external/Lx.mesh", MPI_COMM_WORLD, rank);
  readxml_bcast(ny, docroot, "external/Ly.mesh", MPI_COMM_WORLD, rank);
  readxml_bcast(rseed, docroot, "initialization/Tcvariance/randomseed", MPI_COMM_WORLD, rank);
  readxml_bcast(max_iter, docroot, "solverparameters/Newtonmaxiteration", MPI_COMM_WORLD, rank);
  readxml_bcast(kry_max_iter, docroot, "solverparameters/Krylovmaxiteration", MPI_COMM_WORLD, rank);

  double Ly = 20e-9/LUNIT, tf = 5000e-9/TUNIT, saveperiod = 5.0, alpha = 1e-6, 
         tol_tstep = 0.01, r_t_max = 4.0, r_t_min = 0.2, s_t = 0.9,
         rel_par = 0.9, newton_abs_tol = 1e-6, newton_rel_tol = 1e-3, 
         tramp = 10e-9/TUNIT, sigma = 0.1*double(*TC), Tcvar0 = 2.0*sigma,
         corr_len = 0.2*Ly, TCIMP0 = -20.0/TEMPUNIT, RIMPDIS = 3e-9/LUNIT, DWWIDTH = 10e-9/LUNIT,
         bump = 0.02, T_i = 300.0/TEMPUNIT, eta_i = 0.791296*sqrt(2), mu_i = -0.914352*sqrt(2),
         phi_i = -1000.0, gamma_ei, gamma_hi, curr_i, RVO2_i, Vfrac_i, delV0 = 0.07/VUNIT, delV_i,
         kry_abs_tol = 1e-7, kry_rel_tol = 1e-4;
  readxml_bcast(Ly, docroot, "external/Ly", MPI_COMM_WORLD, rank, LUNIT/1e-9);
  readxml_bcast(tf, docroot, "time/endtime", MPI_COMM_WORLD, rank, TUNIT/1e-9);
  readxml_bcast(delV0, docroot, "external/voltage", MPI_COMM_WORLD, rank, VUNIT);
  readxml_bcast(saveperiod, docroot, "time/saveperiod", MPI_COMM_WORLD, rank);
  readxml_bcast(tramp, docroot, "time/rampt", MPI_COMM_WORLD, rank, TUNIT/1e-9);
  readxml_bcast(T_i, docroot, "initialization/temperature", MPI_COMM_WORLD, rank, TEMPUNIT);
  readxml_bcast(eta_i, docroot, "initialization/SOP", MPI_COMM_WORLD, rank);
  readxml_bcast(mu_i, docroot, "initialization/EOP", MPI_COMM_WORLD, rank);
  readxml_bcast(phi_i, docroot, "initialization/voltage", MPI_COMM_WORLD, rank, VUNIT);
  readxml_bcast(sigma, docroot, "initialization/Tcvariance/sigma", MPI_COMM_WORLD, rank, TEMPUNIT);
  readxml_bcast(corr_len, docroot, "initialization/Tcvariance/correlationlength", MPI_COMM_WORLD, rank, LUNIT/1e-9);
  readxml_bcast(Tcvar0, docroot, "initialization/Tcvariance/mean", MPI_COMM_WORLD, rank, TEMPUNIT);
  readxml_bcast(TCIMP0, docroot, "initialization/Tcvariance/Tcshift", MPI_COMM_WORLD, rank, TEMPUNIT);
  readxml_bcast(RIMPDIS, docroot, "initialization/Tcvariance/radius", MPI_COMM_WORLD, rank, LUNIT/1e-9);
  readxml_bcast(bump, docroot, "initialization/Tcvariance/bump", MPI_COMM_WORLD, rank);
  readxml_bcast(newton_rel_tol, docroot, "solverparameters/Newtonrelativetolerance", MPI_COMM_WORLD, rank);
  readxml_bcast(rel_par, docroot, "solverparameters/Newtonrelaxation", MPI_COMM_WORLD, rank);
  readxml_bcast(tol_tstep, docroot, "solverparameters/timesteptolerance", MPI_COMM_WORLD, rank);
  readxml_bcast(newton_abs_tol, docroot, "solverparameters/Newtonabsolutetolerance", MPI_COMM_WORLD, rank);
  readxml_bcast(kry_abs_tol, docroot, "solverparameters/Krylovabsolutetolerance", MPI_COMM_WORLD, rank);
  readxml_bcast(kry_rel_tol, docroot, "solverparameters/Krylovrelativetolerance", MPI_COMM_WORLD, rank);
  readxml_bcast(alpha, docroot, "solverparameters/Nitschefactor", MPI_COMM_WORLD, rank);

  std::string savemethod("fix"), linear_solver("superlu_dist"), preconditioner("hypre_euclid"),
              Tcvarmethod("random");
  readxml_bcast(savemethod, docroot, "time/savemethod", MPI_COMM_WORLD, rank);
  readxml_bcast(Tcvarmethod, docroot, "initialization/Tcvariance.method", MPI_COMM_WORLD, rank);
  readxml_bcast(linear_solver, docroot, "solverparameters/linearsolver", MPI_COMM_WORLD, rank);
  readxml_bcast(preconditioner, docroot, "solverparameters/preconditioner", MPI_COMM_WORLD, rank);

  enum LogLevel lglvl = INFO;
  readxml_bcast(lglvl, docroot, "solverparameters/loglevel", MPI_COMM_WORLD, rank);
  std::vector<double> t_out;

  //----------------------------------------------------------------
  // Define boundary and initial condition parameter Constant
  //----------------------------------------------------------------
  auto Ts = std::make_shared<Constant>(300.0 / TEMPUNIT);
  readxml_bcast(*Ts, docroot, "external/temperature", MPI_COMM_WORLD, rank, TEMPUNIT);
  auto delV = std::make_shared<Constant>(delV0);                  // V
  auto Resistor = std::make_shared<Constant>(8.0e3 / RUNIT);  // Ohm / [Chosen resistance unit]
  readxml_bcast(*Resistor, docroot, "external/resistor", MPI_COMM_WORLD, rank, RUNIT);
  auto Capacitor = std::make_shared<Constant>(1e-9 / (CUNIT/VUNIT));
  readxml_bcast(*Capacitor, docroot, "external/capacitor", MPI_COMM_WORLD, rank, CUNIT/VUNIT/1e-9);
  auto etas = std::make_shared<Constant>(1.0);
  auto mus = std::make_shared<Constant>(-1.0);
  auto deff = std::make_shared<Constant>(10e-9 / LUNIT);

  if (savemethod == "fix")
  {
    for (double tt = 0.0; tt < tf + 1e-12*tf; tt += saveperiod) t_out.push_back(tt);
  }
  else if (savemethod != "auto")
  {
    std::cout << "User message ===> Warning: input save method not supported, solutions will not be saved" << std::endl;
  }

  gamma_ei = -(double(*CHI) * pow(mu_i, 2)/2 - double(*CHP_IN))/(double(*KB) * T_i);
  gamma_hi = -(double(*CHI) * pow(mu_i, 2)/2 + double(*CHP_IN))/(double(*KB) * T_i);
  RVO2_i = Ly/( double(*ECHARGE) * (double(*NC) * exp(gamma_ei) * double(*MEC) + double(*NV) * exp(gamma_hi) * double(*MHC)) * double(*Lx) * double(*Lz) );
  Vfrac_i = RVO2_i/(RVO2_i + double(*Resistor));
  if (phi_i < -500)
  {
    phi_i = delV0 * Vfrac_i;  // This is default phi_i
    delV_i = delV0;
  }
  else
  {
    delV_i = phi_i / Vfrac_i;
  }
  *delV = delV_i;  // Overloaded operator = for assignment
  curr_i = phi_i / RVO2_i / double(*Lz);

  // For Nitsche's trick dealing with nonstandard boundary conditions
  auto nitsche_eps = std::make_shared<Constant>(alpha*(double(*Lx)/nx));

  MPI_Barrier(MPI_COMM_WORLD);

  //std::cout << "Test node 0" << std::endl;

  //----------------------------------------------------
  // Generate mesh
  //----------------------------------------------------
  // Default SCOTCH is not usable, don't know why, someone online also encountered this and said
  // it has side effects and might be buggy
  parameters["mesh_partitioner"] = "ParMETIS";  
  //std::cout << "Rank " << rank << "'s Lx, Ly, nx, ny: " << double(*Lx) << ", " << Ly << ", " << nx << ", " << ny << std::endl;
  auto mesh = std::make_shared<RectangleMesh>(MPI_COMM_WORLD, Point(0.0, 0.0), Point(double(*Lx), Ly), nx, ny, "crossed");
  //std::cout << "Test node 0.0" << std::endl;
  std::size_t meshdim = mesh->topology().dim();
  // Create mesh functions over the cell facets and initialize to 9
  auto marked_bdrs = std::make_shared<MeshFunction<std::size_t>>(mesh, meshdim - 1, 9);
  // Overloaded operator = for marking all facets as 9
  //*marked_bdrs = 9;

  //std::cout << "Test node 1" << std::endl;

  //--------------------------------------------------------------
  // Create function spaces
  //--------------------------------------------------------------
  auto V = std::make_shared<imt::Form_F::TestSpace>(mesh);
  // Sub FunctionSpace for component of V
  auto V1 = std::make_shared<imt::CoefficientSpace_Tcvar>(mesh);

  //std::cout << "Test node 2" << std::endl;

  //--------------------------------------------------------------
  // Mark boundaries and define standard boundary conditions
  //--------------------------------------------------------------
  auto bias_bdr = std::make_shared<BoundaryY>(0.0);
  auto ground_bdr = std::make_shared<BoundaryY>(Ly);

  bias_bdr->mark(*marked_bdrs, 0);
  ground_bdr->mark(*marked_bdrs, 1);

  auto zero = std::make_shared<Constant>(0.0);
  auto bc_phi_grd = std::make_shared<DirichletBC>(V->sub(4), zero, ground_bdr);
  std::vector<std::shared_ptr<const DirichletBC>> bcs{bc_phi_grd};  // Collect boundary conditions

  //std::cout << "Test node 3" << std::endl;

  //--------------------------------------------------------------
  // Create solution Functions and Coefficients
  //--------------------------------------------------------------
  auto u_n = std::make_shared<Function>(V);
  auto u = std::make_shared<Function>(V);
  auto Tcvar = std::make_shared<TcVariation>(TCIMP0, RIMPDIS, DWWIDTH, double(*Lx));
  // Save Tcvar to file
  File file_Tcvar(figdir+"Tcvar.pvd");
  Function _Tcvar(V1);
  _Tcvar = *Tcvar;   // Overloaded operator = for interpolating Expression to Function
  *(_Tcvar.vector()) *= TEMPUNIT;  // Overloaded operator *= for scaling
  _Tcvar.rename("Tcvar", "Tcvar");
  file_Tcvar << _Tcvar;

  Function ca(V1);  // Charge accumulation
  ca.rename("Charge accumulation", "Charge accumulation");

  // Assign initial conditions and initial guess to solutions
  InitialConditions u_i(eta_i, mu_i, gamma_ei, gamma_hi, phi_i, Ly, T_i, curr_i);
  *u_n = u_i;  // Overloaded operator = for assignment
  *u = u_i;

  //std::cout << "Test node 4" << std::endl;

  //--------------------------------------------------------------------------------
  // Create Forms and attach to them Coefficients, Constants and marked boundaries
  //--------------------------------------------------------------------------------
  auto F = std::make_shared<imt::Form_F>(V);
  auto J = std::make_shared<imt::Form_J>(V, V);
  // integral (s^2 * dx)
  // Seems those formal Forms like ResidualForm etc. don't have their own memory for their coefficients,
  // and the Functions have their pointer passed to the Forms when they are assigned to the Forms, so 
  // the Forms' coefficients will be automatically updated when the external coefficient Functions are
  // updated. But the Functional forms as below seems to have their own memory for their coefficients 
  // and one needs to reassign their coefficients explicitly.
  imt::Form_scalernorm scalernorm(mesh);
  imt::Form_phib0int phi0int(mesh);
  imt::Form_Tint Tint(mesh);
  imt::Form_a a_ca(V1, V1);
  imt::Form_L L_ca(V1);
  Matrix A_ca;
  Vector b_ca;
  assemble(A_ca, a_ca);  

  // Attach marked boundaries for marking boundaries
  F->ds = marked_bdrs;
  J->ds = marked_bdrs;
  phi0int.ds = marked_bdrs;
  L_ca.ds = marked_bdrs;

  // Collect Coefficients and Constants
  std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients_J = 
  {{"dt", dt}, {"Lx", Lx}, {"Lz", Lz}, {"KB", KB}, {"ECHARGE", ECHARGE}, {"TC", TC}, {"T1", T1}, 
   {"AN1", AN1}, {"AN2", AN2}, {"AN3", AN3}, {"T2", T2}, {"AU1", AU1}, {"AU2", AU2}, {"AU3", AU3}, 
   {"GNU1", GNU1}, {"GNU2", GNU2}, {"GNU3", GNU3}, {"CHI", CHI}, {"KN", KN}, {"KAPPAN", KAPPAN}, 
   {"KU", KU}, {"KAPPAU", KAPPAU}, {"NC", NC}, {"NV", NV}, {"MEC", MEC}, {"MEA", MEA}, {"MHC", MHC}, 
   {"MHA", MHA}, {"KEH0", KEH0}, {"PERMITTIVITY", PERMITTIVITY}, {"CPV", CPV}, {"THETA", THETA}, 
   {"HTRAN", HTRAN}, {"CHP_IN", CHP_IN}, {"Resistor", Resistor}, 
   {"Capacitor", Capacitor}, {"deff", deff}, 
   {"nitsche_eps", nitsche_eps}, {"Tcvar", Tcvar}, {"u", u}, {"u_n", u_n}};

  std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients_F = coefficients_J;
  coefficients_F.insert({"etas", etas});
  coefficients_F.insert({"mus", mus});
  coefficients_F.insert({"Ts", Ts});
  coefficients_F.insert({"delV", delV});

  std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients_ca = {{"KB", KB}, 
   {"ECHARGE", ECHARGE}, {"CHI", CHI}, {"NC", NC}, {"NV", NV}, {"MEC", MEC}, {"MEA", MEA}, 
   {"MHC", MHC}, {"MHA", MHA}, {"u", u}};


  // Attach Coefficients to forms (initial conditions and guess included)
  F->set_coefficients(coefficients_F);
  // std::cout << J->coefficient_number("delV") << "/" << J->num_coefficients() << std::endl;  // Test
  J->set_coefficients(coefficients_J);

  L_ca.set_coefficients(coefficients_ca);

  //std::cout << "Test node 5" << std::endl;
  //--------------------------------------------------------------------------------
  // Create nonlinear problem and solver
  //--------------------------------------------------------------------------------
  auto problem = std::make_shared<NonlinearVariationalProblem>(F, u, bcs, J);
  NonlinearVariationalSolver solver(problem);

  // Set solver parameters
  solver.parameters("newton_solver")["relaxation_parameter"] = rel_par;
  solver.parameters("newton_solver")["maximum_iterations"] = max_iter;
  solver.parameters("newton_solver")["absolute_tolerance"] = newton_abs_tol;
  solver.parameters("newton_solver")["relative_tolerance"] = newton_rel_tol;
  solver.parameters("newton_solver")["error_on_nonconvergence"] = false;   // Not stop program if Newton solver does not converge

  solver.parameters("newton_solver")["linear_solver"] = linear_solver; // "superlu_dist", "mumps", "gmres", "bicgstab"
  if (linear_solver != "mumps" && linear_solver != "superlu_dist" && linear_solver != "pastix" 
      && linear_solver != "petsc" && linear_solver != "superlu" && linear_solver != "umfpack")
  {
    solver.parameters("newton_solver")["preconditioner"] = preconditioner;
    solver.parameters("newton_solver")("krylov_solver")["maximum_iterations"] = kry_max_iter;
    solver.parameters("newton_solver")("krylov_solver")["absolute_tolerance"] = kry_abs_tol;
    solver.parameters("newton_solver")("krylov_solver")["relative_tolerance"] = kry_rel_tol;
  }

  // Test
  //double rerror = rel_err(u, u_n, etanorm, munorm, genorm, ghnorm, phinorm, Tnorm);
  //std::cout << "rerror = " << rerror << std::endl;
  //std::cout << "Test node 0" << std::endl;
  //auto usub = std::make_shared<const Function>(*u, 0);
  //std::shared_ptr<Function> s(nullptr);
  //s = split(*u, 0, true);
  ////scalernorm.s = s;
  //scalernorm.set_coefficient("s", s);
  //std::cout << sqrt(assemble(scalernorm)/(double(*Lx)*Ly)) << std::endl;
  //s = split(*u, 1, true);
  //std::cout << sqrt(assemble(scalernorm)/(double(*Lx)*Ly)) << std::endl;
  //std::cout << rel_err(u, u_n, scalernorm) << std::endl;

  File file_eta(figdir+"eta.pvd");
  File file_mu(figdir+"mu.pvd");
  //File file_ge(figdir+"gamma_e.pvd");
  //File file_gh(figdir+"gamma_h.pvd");
  File file_phi(figdir+"phi.pvd");
  File file_T(figdir+"T.pvd");
  File file_ca(figdir+"ca.pvd");
  // Collect data files
  std::vector<File*> files{&file_eta, &file_mu, &file_phi, &file_T, &file_ca};

  //save_sol(files, u, UNITS, 0.0);

  std::ofstream logfile;
  logfile.open("log.txt");
  int pnum = 6, lnum = 15;
  // Set output number format
  std::cout << std::scientific << std::setprecision(pnum);
  logfile << std::scientific << std::setprecision(pnum);
  set_log_level(lglvl);

  //----------------------------------------------------------------------------
  // Start solving the equations using adaptive time stepping
  //----------------------------------------------------------------------------
  double t = 0.0, rerror, dt_factor, mu_norm_av, V_VO2, Tav, R_VO2;
  std::size_t n_step = 0, Tfail = 0, Nfail = 0, otherfail = 0, successive_div = 0;
  std::vector<double>::iterator it_out = t_out.begin();
  std::pair<std::size_t, bool> solvestat;
  auto u2acc = std::make_shared<Function>(V);  // Stores second-order accurate estimate of the solution
  auto dudt_n = std::make_shared<Function>(V);  // Stores du/dt in the previous time step
  auto u_out = std::make_shared<Function>(V);  // Stores linear interpolation between current solution and previous solution
  std::shared_ptr<Function> mu(nullptr);
  mu = split(*u, 1);  // Shallow copy, should automatically update its value after u is updated
  std::shared_ptr<Function> curr(nullptr);

  if (rank == 0)
  {
    logfile << "          #Step            Time       Time step           Tfail           Nfail      Other fail    Av. EOP norm       Av. T (K)  VO2 V drop (V)     VO2 R (Ohm)" << std::endl;
  }

  while (t < tf + 1e-9*tf)
  {
    if (rank == 0)
    {
      std::cout << "------------------------------------------------------------------------" << std::endl;
      std::cout << "User message ===> Start time step refinement at t = " << t 
      << " with dt = " << double(*dt) << " (tol = " << tol_tstep << ")" << std::endl;
    }

    do
    {
      try
      {
        *delV = ramped_delV(t+double(*dt), tramp, delV_i, delV0);  // Overloaded operator = for assignment of the Constant
        //F->delV = delV;  // Update delV of the Form, because Forms have their own memory for their coefficients
        //F->dt = dt;  // Update dt of the Form
        solvestat = solver.solve();
      }
      catch(...)
      {
        otherfail++;
        successive_div++;
        *(u->vector()) = *(u_n->vector());  // Restore initial guess for Newton iteration
        *dt = double(*dt) * r_t_min;   // Decrease dt
        // Mark solve as unsuccessful
        rerror = tol_tstep * 10;  // Set rerror to big value
        dt_factor = 0.1;  // Set dt_factor to small value
        if (rank == 0)
        {
          std::cout << "User message ===> Solver did not converge for other reason! Refined time step dt = " 
          << double(*dt)/r_t_min << " --> " << double(*dt) << " and test again" << std::endl;
        }
        continue;
      }

      if (!solvestat.second)
      {
        Nfail++;
        successive_div++;
        *(u->vector()) = *(u_n->vector());  // Restore initial guess for Newton iteration
        *dt = double(*dt) * r_t_min;   // Decrease dt
        // Mark solve as unsuccessful
        rerror = tol_tstep * 10;  // Set rerror to big value
        dt_factor = 0.1;  // Set dt_factor to small value
        if (rank == 0)
        {
          std::cout << "User message ===> Newton solver did not converge! Refined time step dt = " 
          << double(*dt)/r_t_min << " --> " << double(*dt) << " and test again" << std::endl;
        }
        continue;
      }

      successive_div = 0;  // Set back to zero because solve is successful

      // Compute the second-order accurate estimate of the solution; the current solution u is first-order accurate (backward Euler)
      if (n_step > 0)
      {
        // u2acc = u_n + (dudt_n + (u - u_n)/dt)*dt/2
        *(u2acc->vector()) = *(u_n->vector()); 
        u2acc->vector()->axpy(double(*dt)/2, *(dudt_n->vector()));
        u2acc->vector()->axpy(0.5, *(u->vector()));
        u2acc->vector()->axpy(-0.5, *(u_n->vector()));
      }
      else  // Mark solve as successful but dt_factor as small
      {  
        rerror = tol_tstep * 0.1;  // Set rerror to small value
        dt_factor = 0.1;  // Set dt_factor to small value
        break;
      }

      // Calculate the backward Euler time integration error and time step changing factor
      rerror = rel_err(u2acc, u, scalernorm);
      dt_factor = std::min(std::max(s_t*sqrt(tol_tstep/std::max(rerror, 1e-10)), r_t_min), r_t_max);

      // Adjust things for adaptive time stepping
      if (rerror > tol_tstep)
      {
        Tfail++;
        *(u->vector()) = *(u_n->vector());  // Restore initial guess for Newton iteration
        *dt = double(*dt) * dt_factor;
        if (rank == 0)
        {
          std::cout << "User message ===> Time stepping error " << rerror 
          << " is too big, refined time step dt = " << double(*dt)/dt_factor 
          << " --> " << double(*dt) << " and test again" << std::endl;
        }
      }
    } while (rerror > tol_tstep && successive_div <= max_div);

    if (successive_div > max_div) break;

    // Time step refinement is done, now save solution, update solution, and output
    n_step++;
    if (savemethod == "auto")
    {
      if (n_step % lround(saveperiod) == 1) {
        L_ca.u = u;
        assemble(b_ca, L_ca);
        solve(A_ca, *(ca.vector()), b_ca);  // Cast charge accumulation Function
        save_sol(files, u, ca, UNITS, t+double(*dt));
      }
    }
    else if (savemethod == "fix")
    {
      while (it_out < t_out.end() && *it_out <= t+double(*dt))
      {
        // Compute linear interpolation between current solution and previous solution
        // u_out = u_n + (u-u_n)*(t_out[n_out]-t)/dt
        *(u_out->vector()) = *(u_n->vector());
        u_out->vector()->axpy((*it_out-t)/double(*dt), *(u->vector()));
        u_out->vector()->axpy(-(*it_out-t)/double(*dt), *(u_n->vector()));
        L_ca.u = u_out;
        assemble(b_ca, L_ca);
        solve(A_ca, *(ca.vector()), b_ca);  // Cast charge accumulation Function
        save_sol(files, u_out, ca, UNITS, *it_out);
        if (rank == 0)
        {
          std::cout << "User message ===> Out: " << it_out-t_out.begin() << ", t = " << *it_out << std::endl;

          logfile << " - " << std::setw(lnum-3) << it_out-t_out.begin() << " " << std::setw(lnum) 
          << *it_out << " - out" << std::endl;
        }
        it_out++;
      }
    }

    scalernorm.s = mu;
    mu_norm_av = sqrt(assemble(scalernorm)/(double(*Lx)*Ly));  // Calculate norm average of mu
    phi0int.u = u;
    V_VO2 = assemble(phi0int)/double(*Lx)*VUNIT;  // Calculate voltage drop across VO2
    Tint.u = u;
    Tav = assemble(Tint)/(double(*Lx)*Ly)*TEMPUNIT;  // Calculate average temperature across VO2
    curr = split(*u, 6, true);  // Deep copy, used to access the subfunciton's dof
    R_VO2 = V_VO2 / (double(*Lz) * (*(curr->vector()))[0] * CUNIT/TUNIT);  // # Calculate VO2 resistance
    if (rank == 0)
    {
      std::cout << "User message ===> Completed refinement: dt = " << double(*dt) 
      << ", t = " << t << " --> " << t+double(*dt) << std::endl;

      logfile << std::setw(lnum) << n_step << " " << std::setw(lnum) << t+double(*dt) << " " << std::setw(lnum) 
      << double(*dt) << " " << std::setw(lnum) << Tfail << " " << std::setw(lnum) << Nfail << " " 
      << std::setw(lnum) << otherfail << " " << std::setw(lnum) << mu_norm_av << " " << std::setw(lnum) << Tav 
      << " " << std::setw(lnum) << V_VO2 << " " << std::setw(lnum) << R_VO2 << std::endl;
    }

    // Update dudt_n = (u - u_n)/dt
    *(dudt_n->vector()) = *(u->vector()); *(dudt_n->vector()) /= double(*dt);
    dudt_n->vector()->axpy(-1/double(*dt), *(u_n->vector()));
    *(u_n->vector()) = *(u->vector());  // Update the previous solution
    t += double(*dt);

    if (dt_factor > 1)
    {
      *dt = double(*dt) * dt_factor;
      if (rank == 0)
      {
        std::cout << "User message ===> dt is increased: dt = " << double(*dt)/dt_factor 
        << " --> " << double(*dt) << std::endl;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  std::clock_t c_end = std::clock();  // Record the calculation end time
  double time_elapsed = (c_end-c_start)/CLOCKS_PER_SEC;

  if (successive_div > max_div)
  {
    if (rank == 0)
    {
      std::cout << "User message ===> !!! Solving process diverged, stop the process !!! Computation time: " 
      << time_elapsed << " s." << std::endl;

      logfile << "!!! Solving process diverged, stop the process !!! Computation time: " 
      << std::setw(lnum) << time_elapsed << " s." << std::endl;
    }
  }
  else
  {
    if (rank == 0)
    {
      std::cout << "User message ===> Finished computation, computation time: " 
      << time_elapsed << " s." << std::endl;
      
      logfile << "Finished computation, computation time: " 
      << std::setw(lnum) << time_elapsed << " s." << std::endl;
    }
  }

  logfile.close();

  return 0;
}
