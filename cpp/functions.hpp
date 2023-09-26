#pragma once

#include <cmath>
#include <math.h> 
#include <complex>
#include <vector>
#include <stdint.h>
#include <cstdint>

typedef double data_t;
#define EPS 1e-10


// Build the necessary matrices and their compound matrices (for P-SV plane waves)
// Matrices are stored in row major order

void Q_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t t); // = M.P
void invQ_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t g); // = P-1.M-1
void field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z); // = Q.E
void inv_field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z); // = E-1.Q-1
void propagator_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z); // = Q.E.Q-1
void compound_Q_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t t);
void compound_invQ_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t g);
void compound_field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z);
void compound_field_matrix_c(std::complex<data_t> * M, const data_t f, const std::complex<data_t> k, const data_t a, const data_t b, const data_t rho, const data_t z);
void compound_inv_field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z);
void compound_inv_field_matrix_c(std::complex<data_t> * M, const data_t f, const std::complex<data_t> k, const data_t a, const data_t b, const data_t rho, const data_t z);
void compound_propagator_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z);
void compound_propagator_matrix_c(std::complex<data_t> * M, const data_t f, const std::complex<data_t> k, const data_t a, const data_t b, const data_t rho, const data_t z);


// Compute the characteristic function for a given configuration of guided waves
// v=(A+,A-,B+,B-)' is the amplitudes vector ; + for downgoing, - for upgoing, A for P-waves, B for S-waves
// w=(ux,uz,sxz,szz)' is the state variables vector; u is the displacement, s is the stress
// w = M.v  M is a 4 x 4 matrix, product of field and propagator matrices

// P-SV Lamb waves
// M = T; 0 = T(3:4,1:2).(ux,uz)'
std::complex<data_t> free_solid_free(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);
void free_solid_free_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);

// P-SV Rigid walls
// M = T; 0 = T(1:2,3:4).(sxz,szz)'
std::complex<data_t> rigid_solid_rigid(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);
void rigid_solid_rigid_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);

// P-SV Surface waves
// M = T.D; 0 = T(3:4,:).D(:,1&3)
std::complex<data_t> free_solid_halfspace(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);
void free_solid_halfspace_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);

// P-SV Embedded guided waves
// M = D^-1.T.D; 0 = D^-1(2&4,:).T(:,:).D(:,2&4)
std::complex<data_t> halfspace_solid_halfspace(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d); 
void halfspace_solid_halfspace_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);

// P-SV Embedded guided waves - leaky modes
std::complex<data_t> halfspace_solid_halfspace_leaky(const data_t f, const std::complex<data_t> k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d); 

// P-SV Surface waves - leaky modes
std::complex<data_t> free_solid_halfspace_leaky(const data_t f, const std::complex<data_t> k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d);

// customized square root for complex numbers
std::complex<data_t> sqrt_mb(std::complex<data_t> x);

// print function for QC
void print(std::complex<data_t> * M, int m, int n);