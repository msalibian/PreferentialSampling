#include <TMB.hpp>
#include<cmath>
//#include <tmb_core.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_VECTOR(Y1);        // Sample location longitudes. Must have same length as Y
  DATA_VECTOR(Y2);        // Sample location latitudes. Must have same length as Y
  DATA_VECTOR(Y);         // Observations
  DATA_VECTOR(Ind);         // Used to create mean vector
  DATA_IVECTOR(pointer);   // points to rows of S which are closest to sampling locations
  DATA_IVECTOR(meshidxloc);

  PARAMETER_VECTOR(S); // latent field to be integrated out

  DATA_STRUCT(spde,spde_t);
  // Parameters
  PARAMETER(mu);
  PARAMETER(log_phi);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  PARAMETER(beta);
  // Do some internal reparameterisation
  Type kappa = exp(log_kappa);
  Type phi = exp(log_phi);
  Type tau = exp(log_tau);

  Type sigma = 1 / sqrt(4 * M_PI * exp(-2*log_phi) * exp(2*log_kappa));

  int i;
  Type ans=0;  // ans will be the resulting likelihood.

  Type intSum1=0;
  // First chunk of likelihood
  // initialise precision matrix and mean vector
  // [S]
  SparseMatrix<Type> Q = Q_spde(spde,(1/phi));

  vector<Type> muvec(S.size());
  muvec=Ind*mu;
  ans += SCALE(GMRF(Q), (1/kappa))(S - muvec);      // Negative log likelihood

  // Second chunk of likelihood
  // [Y|S, X]

  Type intSum2=0;
  for (i=0;i<Y.size();i++)
  {
    ans -= dnorm(Y(i), S(pointer(i)), tau, 1);
    // calculate part of third chunk below
    intSum2 += S(pointer(i));
  }
  // Third chunk of likelihood
  // [X|S]
  // calculate denominator  int exp(S(u)) du
  vector<Type> preSum(meshidxloc.size());
  preSum=beta/dexp(S(meshidxloc), Type(beta), 0);
  intSum1 += sum(preSum);
  // multiply integral sum by the area of each lattice square
  Type logInt=0;
  logInt += log(intSum1*(0.03225806*0.03225806));

  ans -= (beta*intSum2) - (pointer.size()*logInt);

  REPORT(sigma);
  return ans;
  }
