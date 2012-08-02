#include <dolfin.h>   //Include the FEniCS code
#include "utils.h"
#include "Poisson.h"
#include "L2error.h"
#include "Velocity.h"
#include "OmegaAD.h"
#include "ThetaAD.h"

using namespace dolfin;

//Initial Temperature Profile (uniform stratification)
class Tini : public Expression
{
public:

  Tini(const double& N2, const double& grav, const double& alphat, const double& Tmean)
      : Expression(), N2(N2), grav(grav), alphat(alphat), Tmean(Tmean) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    const double Tz = N2 / (grav * alphat);
    values[0] = Tz*(1000.0 * x[1] - 500) + Tmean;
  }

private:

   const double& N2;
   const double& grav;
   const double& alphat;
   const double& Tmean;
};

//Define Diffusivity Ratio (scaled-Vertical Diffusivity)
class DiffRatio : public Expression
{
public:

 DiffRatio(const double& kappaL, const double& kappaR, const double& efold, const double& yRise, const double& mRatio, const double& kappah, const double& delPeak, const double& delBdd)
		: Expression(), kappaL(kappaL), kappaR(kappaR), efold(efold), yRise(yRise), mRatio(mRatio), kappah(kappah), delPeak(delPeak), delBdd(delBdd) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {

	double kMax;
   double kappav;
   double D; 

   //Catch the const kMax case before calling function 
	if (kappaL == kappaR)
	{	kMax = kappaL + 0.*x[0]; }
	else
	{	kMax = kappaMax(x[0],kappaL,kappaR,delBdd,delPeak); }	

   D = std::max(depth(x[0],yRise,mRatio),0.);

   kappav = kMax*exp(-1000*(x[1]-D)/efold);
   
   //Ratio is kappav/kappah 
   values[0] = kappav/kappah;

  }

private:

  const double& kappaL;
  const double& kappaR; 
  const double& efold;
  const double& yRise;
  const double& mRatio;
  const double& kappah;
  const double& delPeak;
  const double& delBdd;
};

//Define Subdomain map for periodic boundary condition
class PeriodicBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] > -DOLFIN_EPS && x[0] < DOLFIN_EPS && on_boundary;
  }

  void map(const Array<double>& x, Array<double>& y) const
  {
    y[0] = x[0] - 2.0;
    y[1] = x[1];
  }
};

// ********* MAIN ***********
int main(int argc, char *argv[])
{
  init(argc, argv);
  parameters["allow_extrapolation"] = true;

  logging(false);


//cout << "No text?" << endl;  
//Read mesh
  std::string inMesh, inSubD, outMesh;
  inMesh.reserve(PETSC_MAX_PATH_LEN);
  GetFilename("-mesh",inMesh);
  inSubD.reserve(PETSC_MAX_PATH_LEN);
  GetFilename("-subd",inSubD);
//  outMesh.reserve(PETSC_MAX_PATH_LEN);
//  GetFilename("-outM",outMesh);
  Mesh mesh(inMesh);
//  Mesh oMesh(outMesh);
  MeshFunction<unsigned int> sub_domains(mesh, inSubD);
cout << "Text printing works" << endl;

//Define function space
  OmegaAD::FunctionSpace V(mesh);

//Base parameters 
  double   vish, visv, difh, H_val, L, g, alphat, nbc_val;
  double   vh2, H3;
  vish    = 1E-3;
   vh2    = vish*vish;
  visv    = 1E-3;
  difh    = 1E-3;
  H_val   = 1000;
   H3	  = H_val*H_val*H_val;
  L       = 100000;
  g       = 9.8;
  alphat  = 2E-4;
  nbc_val = 0.0;

  PetscReal yRise = 400;
  PetscOptionsGetReal(PETSC_NULL,"-yRise",&yRise,PETSC_NULL);
  PetscReal kvR = 0.1;
  PetscOptionsGetReal(PETSC_NULL,"-kvR",&kvR,PETSC_NULL);
  PetscReal kvL = kvR;
  PetscOptionsGetReal(PETSC_NULL,"-kvL",&kvL,PETSC_NULL);
  PetscReal mRatio = 1.0;
  PetscOptionsGetReal(PETSC_NULL,"-mRatio",&mRatio,PETSC_NULL);
  PetscReal delPeak = 6000;
  PetscOptionsGetReal(PETSC_NULL,"-delPeak",&delPeak,PETSC_NULL);
  PetscReal delBdd = 10000;
  PetscOptionsGetReal(PETSC_NULL,"-delBdd",&delBdd,PETSC_NULL);
cout << "Base params set" << endl;

  PetscReal N2freq = 1E-3;
  PetscOptionsGetReal(PETSC_NULL,"-N2freq",&N2freq,PETSC_NULL);
  PetscReal eFold = 250;
  PetscOptionsGetReal(PETSC_NULL,"-eFold",&eFold,PETSC_NULL);

//Nondimensionals
  Constant    Ar(H_val/L);
  Constant   iPr(difh/vish);
  Constant rvisc(visv/vish);
  Constant    Ra(g*alphat*H3/vh2);
  Constant   nbc(nbc_val);
  Constant     H(H_val);
 
  DiffRatio rk(kvL, kvR, eFold, yRise, mRatio, difh, delPeak/L, delBdd/L);
  Function rdiff(V);
    rdiff.interpolate(rk);
cout << "Nondims set" << endl;

  printf("Aspect ratio	= %g\n", H/L);
  printf("Prandtl num	= %g\n", vish/difh);
  printf("R_visc	= %g\n", visv/vish);
  printf("Rayleigh num  = %g\n", g*alphat*H3/vh2);
  
//Timestepping parameters
  PetscReal        dt=600;
  PetscReal      days=1.;
  PetscInt    dt_plot=43200;
  PetscOptionsGetReal(PETSC_NULL,"-dt",&dt,PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL,"-t_days",&days,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-dt_plot",&dt_plot,PETSC_NULL);

  PetscReal    dt_max=dt;
  PetscOptionsGetReal(PETSC_NULL,"-dt_max",&dt_max,PETSC_NULL);
  Constant       k(dt*1E-11);

// Adaptive parameters
  double tol = 1E-2;
  double normerr = 0.;
 
  char namestr[100];
  std::string filename; 
 
  dolfin::uint n=0;
  double       t=0;
  double       T=86400*days;	//Times to hour-halfday-day-#

  
  //Initiate Functions
  Tini   T_ini(N2freq, 9.81, 2E-4, 3.0); 	// Stratified temp(N2, g, alphat, Tmean)
  Function  th1(V);
  Function  th0(V);
   th0.interpolate(T_ini);
  Function  om1(V);
  Function  om0(V);
   om0.vector().zero();
  Function psi0(V);
  Function psi1(V);
   psi1.vector().zero();
  Function psin(V);
  Function z0(V);
   z0.vector().zero();

  //----BC----
  Constant bd_0(0.0);
  DirichletBC bct0(V, bd_0, sub_domains, 0);
  DirichletBC bcr0(V, bd_0, sub_domains, 1);
  DirichletBC bcb0(V, bd_0, sub_domains, 2);
  DirichletBC bcl0(V, bd_0, sub_domains, 3);
  
  DirichletBC bct_th(V, T_ini, sub_domains, 0);
  
  PeriodicBoundary periodic_bdd;
  PeriodicBC pbc(V, periodic_bdd);
 
  std::vector<const BoundaryCondition*> bcs_psi;
  std::vector<const BoundaryCondition*> bcs_th;
  std::vector<const BoundaryCondition*> bcs_om;

  bcs_psi.push_back(&bct0);
  bcs_psi.push_back(&bcb0);
  bcs_psi.push_back(&bcl0);
  bcs_psi.push_back(&bcr0);
  //bcs_psi.push_back(&pbc);
  
  bcs_om.push_back(&bcb0);
  bcs_om.push_back(&bct0);
  bcs_om.push_back(&bcl0);
  bcs_om.push_back(&bcr0);
  //bcs_om.push_back(&pbc);
  
  //bcs_th.push_back(&bct_th);
  //bcs_th.push_back(&pbc); 

  //------ Theta -----
  ThetaAD::BilinearForm a_th(V, V);
  ThetaAD::LinearForm L_th(V);
  a_th.k   = k;	    a_th.Ar = Ar; 
  a_th.iPr = iPr; a_th.r  = rdiff;
  a_th.psi1= psin;
  L_th.k   = k;     L_th.Ar = Ar;    L_th.nbc = nbc;
  L_th.iPr = iPr; L_th.r  = rdiff;
  L_th.psi0= psi0;  L_th.th0= th0;   
  VariationalProblem prob_th(a_th, L_th, bcs_th);
cout << "Setup Temp";

  //------ Omega -----
  OmegaAD::BilinearForm a_om(V, V);
  OmegaAD::LinearForm L_om(V);
  a_om.k   = k;     a_om.Ar = Ar;
  a_om.r  = rvisc;
  a_om.psi1= psin;
  L_om.k   = k;	    L_om.Ar = Ar;    L_om.nbc = nbc;
  L_om.r  = rvisc; L_om.Ra = Ra;
  L_om.psi0= psi0;  L_om.th0= th0;   L_om.om0= om0;   L_om.th1= th1; 
  VariationalProblem prob_om(a_om, L_om, bcs_om);
cout << " and Vort";  

  //------ Psi -------
  Poisson::BilinearForm a_psi(V, V);
  Poisson::LinearForm L_psi(V);
  a_psi.Ar = Ar;
  L_psi.om1 = om1;  L_psi.nbc = nbc;
  VariationalProblem prob_psi(a_psi, L_psi, bcs_psi);
cout << " and Streamfunction"; 

  //----- L2norm -----
  L2error::Functional L2(mesh);
  L2.u = psin; L2.u_true = psi1;
cout << "and norm\n";

  //------ Output ----
  Mesh oMesh(genMesh(yRise, mRatio));
  File smeshFile("out/squareMesh.xml");
  File smeshFilepvd("out/squareMesh.pvd");
  smeshFile << oMesh;
  smeshFilepvd << oMesh;
  
	Poisson::FunctionSpace Vout(oMesh);
	Function psiout(Vout);
	Function thout(Vout);
//	Function omout(Vout);
	psiout.vector().zero();
	thout.vector().zero();
//	omout.vector().zero();

  File file_th("out/temp.pvd", "compressed");
  File file_om("out/vorticity.pvd", "compressed");
  File file_psi("out/psi.pvd", "compressed");
  File file_kv("out/kappav.pvd", "compressed");
   file_th  << th0;
   file_om  << om0;
   file_psi << psi0;
   file_kv  << rdiff;



//***********************************************
//*****   Loop
//***********************************************

cout << "Entering loop" << endl;
dolfin_debug("Entering Loop");
  while (t <= T && n < 20)
  {
     n = 0; psin = psi0;
     while(n < 2 || (normerr > tol && n < 21))
     { 
       psin = psi1;
       prob_th.solve(th1);
       prob_om.solve(om1);
       prob_psi.solve(psi1);
       normerr = assemble(L2);	//  ||psi1 - psin||_2
       n += 1;
     };

     printf(">> Convergence difference of %f after %i steps, at time %g\n",normerr,n,t);

     //Dump when its taking more than 10 iterations
     if (fmod(t,dt_plot) == 0 || n > 10) {
        printf(">> Saving results\n\n");
			file_th  << th1;
			file_om  << om1;
			file_psi << psi1;

      sprintf(namestr,"out/psi_t%.0f.xml",t);
      filename.assign(namestr);
 		File psiFile(filename);
			psiout.interpolate(psi1);
			psiFile << psiout.vector();
		sprintf(namestr,"out/theta_t%.0f.xml",t);
		filename.assign(namestr);
		File thFile(filename);
			thout.interpolate(th1);
			thFile << thout.vector();

     };

     // Update functions
     th0 = th1;
     om0 = om1;
     psi0= psi1;
     t += dt;
  }

}
