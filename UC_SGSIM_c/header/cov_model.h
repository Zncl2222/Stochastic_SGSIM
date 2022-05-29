
# ifndef COV_MODEL_H
# define COV_MODEL_H

void Cov_model(double *x, double* cov,int n_dim, double a, double C0);

void Cov_model2d(double **x,double* cov,int n_dim, double a, double C0);

# endif