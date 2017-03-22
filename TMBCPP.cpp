#include <TMB.hpp>
#include<cmath>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y1);        // Sample location longitudes. Must have same length as Y
  DATA_VECTOR(Y2);        // Sample location latitudes. Must have same length as Y
  DATA_VECTOR(Y);         // Observations
  DATA_MATRIX(distMat);   // Distances between observations
  DATA_IVECTOR(pointer);   // points to rows of S which are closest to sampling locations

  PARAMETER_VECTOR(S);

  PARAMETER(mu);
  PARAMETER(phi);
  PARAMETER(sigma);
  PARAMETER(tau);
  PARAMETER(beta);

  using namespace density;
  int i,j;
  Type ans=0;  // ans will be the resulting likelihood
  Type intSum1=0;
  // First chunk of likelihood
  // initialise covariance matrix and mean vector
  matrix<Type> cov(S.size(),S.size());
  vector<Type> muvec(S.size());
  for (i=0;i<S.size();i++)
  {
    muvec(i) = mu;
    cov(i,i)=exp(sigma)*exp(sigma);
    intSum1 += exp(beta*S(i));
    for ( j=0;j<i;j++)
    {
      cov(i,j)=matern(distMat(i,j), exp(phi), Type(1))*(exp(sigma)*exp(sigma));
      cov(j,i)=cov(i,j);
    }
  }
  density::MVNORM_t<Type> neg_log_density(cov);
  // [S]
  ans += neg_log_density(S-muvec);
  // Second chunk of likelihood
  // initialise conditional covariance and mean
  matrix<Type> covcond(Y.size(),Y.size());
  vector<Type> muveccond(Y.size());
  for (i=0;i<Y.size();i++)
  {
    muveccond(i) = S(pointer(i));
    covcond(i,i)=exp(tau)*exp(tau);
    for ( j=0;j<i;j++)
    {
      covcond(i,j)=0;
      covcond(j,i)=covcond(i,j);
    }
  }
  density::MVNORM_t<Type> neg_log_density_cond(covcond);
  // [Y|S, X]
  ans += neg_log_density_cond(Y-muveccond);
  // Third chunk of likelihood
  // initialise the sampling likelihood (Poisson Process)
  Type logInt=0;
  logInt += log(intSum1*(0.05*0.05));

  Type intSum2=0;
  for(int i=0;i<pointer.size();i++){
    intSum2 += S(pointer(i));
  }

  Type logLik=0;
  logLik += (beta*intSum2) - (pointer.size()*logInt);
  // [X|S]
  ans -= logLik;
  return ans;
  }
