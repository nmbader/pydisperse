#include "functions.hpp"
#include <iostream>
#include <memory>
#include <omp.h>

int sign(data_t x){
    if (x < 0) return -1;
    else return 1;
}

void Q_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t t){
    M[0]=1;
    M[1]=1;
    M[2]=-s;
    M[3]=s;
    M[4]=r;
    M[5]=-r;
    M[6]=1;
    M[7]=1;
    M[8]=-2*mu*r;
    M[9]=2*mu*r;
    M[10]=-mu*t;
    M[11]=-mu*t;
    M[12]=mu*t;
    M[13]=mu*t;
    M[14]=-2*mu*s;
    M[15]=2*mu*s;
}

void invQ_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t g){
    M[0]=g;
    M[1]=(1-2*g)/(2.0*r);
    M[2]=-g/(2*mu*r);
    M[3]=-g/(2*mu);
    M[4]=g;
    M[5]=-M[1];
    M[6]=-M[2];
    M[7]=M[3];
    M[8]=(2*g-1)/(2.0*s);
    M[9]=g;
    M[10]=-M[3];
    M[11]=-g/(2*mu*s);
    M[12]=-M[8];
    M[13]=g;
    M[14]=-M[3];
    M[15]=-M[11];
}

void field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    data_t c=2*M_PI*f/k;
    std::complex<data_t> r=sqrt((std::complex<data_t>)((c/a)*(c/a)-1));
    std::complex<data_t> s=sqrt((std::complex<data_t>)((c/b)*(c/b)-1));
    data_t t=2-(c/b)*(c/b);    
    std::complex<data_t> er=exp(i*k*r*z);
    std::complex<data_t> es=exp(i*k*s*z);
    M[0] = er;
    M[1] = 1.0/er;
    M[2] =-s*es;
    M[3] =s/es;
    M[4] =r*er;
    M[5] =-r/er;
    M[6] =es;
    M[7] =1.0/es;
    M[8] =-2*mu*r*er;
    M[9] =2*mu*r/er;
    M[10]=-mu*t*es;
    M[11]=-mu*t/es;
    M[12]=mu*t*er;
    M[13]=mu*t/er;
    M[14]=-2*mu*s*es;
    M[15]=2*mu*s/es;
}

void inv_field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    data_t c=2*M_PI*f/k;
    std::complex<data_t> r=sqrt((std::complex<data_t>)((c/a)*(c/a)-1));
    std::complex<data_t> s=sqrt((std::complex<data_t>)((c/b)*(c/b)-1));
    data_t g=(b/c)*(b/c);
    std::complex<data_t> er=exp(i*k*r*z);
    std::complex<data_t> es=exp(i*k*s*z);
    M[0] =g/er;
    M[1] =(1-2*g)/(2.0*r*er);
    M[2] =-g/(2*mu*r*er);
    M[3] =-g/(2*mu*er);
    M[4] =g*er;
    M[5] =(2*g-1)/(2.0*r)*er;
    M[6] =g/(2*mu*r)*er;
    M[7] =-g/(2*mu)*er;
    M[8] =(2*g-1)/(2.0*s*es);
    M[9] =g/es;
    M[10]=g/(2*mu*es);
    M[11]=-g/(2*mu*s*es);
    M[12]=-(2*g-1)/(2.0*s)*es;
    M[13]=g*es;
    M[14]=g/(2*mu)*es;
    M[15]=g/(2*mu*s)*es;
}

void propagator_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    data_t c=2*M_PI*f/k;
    std::complex<data_t> r=sqrt((std::complex<data_t>)((c/a)*(c/a)-1));
    std::complex<data_t> s=sqrt((std::complex<data_t>)((c/b)*(c/b)-1));
    data_t g=(b/c)*(b/c);
    data_t t=2-(c/b)*(c/b);
    std::complex<data_t> Ca=cos(k*r*z);
    std::complex<data_t> Cb=cos(k*s*z);
    std::complex<data_t> Sa=sin(k*r*z);
    std::complex<data_t> Sb=sin(k*s*z);
    M[0] =g*(2.0*Ca-t*Cb);
    M[1] =-i*g*(t/r*Sa+2.0*s*Sb);
    M[2] =-i*g/mu*(1.0/r*Sa+s*Sb);
    M[3] =-g/mu*(Ca-Cb);
    M[4] =i*g*(2.0*r*Sa+t/s*Sb);
    M[5] =-g*(t*Ca-2.0*Cb);
    M[6] =M[3];
    M[7] =-i*g/mu*(r*Sa+1.0/s*Sb);
    M[8] =-i*g*mu*(4.0*r*Sa+t*t/s*Sb);
    M[9] =2*g*mu*t*(Ca-Cb);
    M[10]=M[0];
    M[11]=M[4];
    M[12]=M[9];
    M[13]=-i*g*mu*(t*t/r*Sa+4.0*s*Sb);
    M[14]=M[1];
    M[15]=M[5];
}

void compound_Q_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t t){
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    data_t mt=mu*t;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    M[0] =-2.0*r;   M[1] =1.0+rs;          M[2] =1.0-rs;           M[3] =M[2];        M[4] =M[1];  M[5] =-2.0*s;
    M[6] =4.0*mr;   M[7] =-mt-2.0*mrs;     M[8] =-mt+2.0*mrs;      M[9] =M[8];        M[10]=M[7];  M[11]=2*mt*s;
    M[12]=0;        M[13]=-2.0*ms+mt*s;    M[14]=-M[13];           M[15]=M[13];       M[16]=-M[13];M[17]=0;
    M[18]=0;        M[19]=-mt*r+2.0*mr;    M[20]=M[19];            M[21]=mt*r-2.0*mr; M[22]=M[21]; M[23]=0;
    M[24]=2*mt*r;   M[25]=-2.0*mrs-mt;     M[26]=2.0*mrs-mt;       M[27]=M[26];       M[28]=M[25]; M[29]=4.0*ms;
    M[30]=-4*mt*mr; M[31]=4.0*mr*ms+mt*mt; M[32]=-4.0*mr*ms+mt*mt; M[33]=M[32];       M[34]=M[31]; M[35]=-4*mt*ms;
}

void compound_invQ_matrix(std::complex<data_t> * M, const data_t mu, const std::complex<data_t> r, const std::complex<data_t> s, const data_t g){
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    M[0] =g*(2*g-1)/r;                      M[1] =g*g/mr;                           M[2] =0;            M[3] =0;            M[4] =M[0]/(2*mu);   M[5] =M[1]/(2*mu);
    M[6] =g*g+(2.0*g-1)*(2.0*g-1)/(4.0*rs); M[7] =g*g/(2*mu)+g*(2*g-1)/(4.0*mrs);   M[8] =-g/(4.0*ms);  M[9] =g/(4.0*mr);   M[10]=M[7];          M[11]=(g/(2*mu))*(g/(2*mu))*(1.0+1.0/rs);
    M[12]=2*g*g-M[6];                       M[13]=g*g/mu-M[7];                      M[14]=-M[8];        M[15]=M[9];         M[16]=M[13];         M[17]=(g/(2*mu))*(g/(2*mu))*(1.0-1.0/rs);
    M[18]=M[12];                            M[19]=M[13];                            M[20]=M[8];         M[21]=-M[9];        M[22]=M[13];         M[23]=M[17];
    M[24]=M[6];                             M[25]=M[7];                             M[26]=-M[8];        M[27]=-M[9];        M[28]=M[7];          M[29]=M[11];
    M[30]=g*(2*g-1)/s;                      M[31]=M[30]/(2*mu);                     M[32]=0;            M[33]=0;            M[34]=g*g/ms;        M[35]=M[34]/(2*mu);
}

void compound_field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    data_t c=2*M_PI*f/k;
    if ((c == a) || (c == b)) c += EPS;
    std::complex<data_t> r=sqrt((std::complex<data_t>)((c/a)*(c/a)-1));
    std::complex<data_t> s=sqrt((std::complex<data_t>)((c/b)*(c/b)-1));
    data_t t=2-(c/b)*(c/b);    
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    data_t mt=mu*t;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    std::complex<data_t> e2=exp(i*k*(r+s)*z);
    std::complex<data_t> e3=exp(i*k*(r-s)*z);
    std::complex<data_t> e4=exp(-i*k*(r-s)*z);
    std::complex<data_t> e5=exp(-i*k*(r+s)*z);
    M[0] =-2.0*r;   M[1] =1.0+rs;            M[2] =1.0-rs;              M[3] =M[2];         M[4] =M[1];     M[5] =-2.0*s;
    M[6] =4.0*mr;   M[7] =-mt-2.0*mrs;       M[8] =-mt+2.0*mrs;         M[9] =M[8];         M[10]=M[7];     M[11]=2*mt*s;
    M[12]=0;        M[13]=-2.0*ms+mt*s;      M[14]=2.0*ms-mt*s;         M[15]=M[13];        M[16]=M[14];    M[17]=0;
    M[18]=0;        M[19]=-mt*r+2.0*mr;      M[20]=M[19];               M[21]=mt*r-2.0*mr;  M[22]=M[21];    M[23]=0;
    M[24]=2*mt*r;   M[25]=-2.0*mrs-mt;       M[26]=2.0*mrs-mt;          M[27]=M[26];        M[28]=M[25];    M[29]=4.0*ms;
    M[30]=-4*mt*mr; M[31]=4.0*mr*ms+mt*mt;   M[32]=-4.0*mr*ms+mt*mt;    M[33]=M[32];        M[34]=M[31];    M[35]=-4*mt*ms;
    for (int j=0; j<6; j++){
        M[j*6+1] *= e2;
        M[j*6+2] *= e3;
        M[j*6+3] *= e4;
        M[j*6+4] *= e5;
    }
}

void compound_field_matrix_c(std::complex<data_t> * M, const data_t f, const std::complex<data_t> k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    std::complex<data_t> c=2*M_PI*f/k;
    if ((c.real() == a) || (c.real() == b)) c += EPS;
    std::complex<data_t> r=sqrt_mb((c/a)*(c/a)-1.0);
    std::complex<data_t> s=sqrt_mb((c/b)*(c/b)-1.0);
    //r *= (data_t)sign(r.imag());
    //r *= (data_t)sign(k.real()*r.real()-k.imag()*r.imag());
    //s *= (data_t)sign(s.imag());
    //s *= (data_t)sign(k.real()*s.real()-k.imag()*s.imag());
    std::complex<data_t> t=-(c/b)*(c/b)+2.0;    
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    std::complex<data_t> mt=mu*t;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    std::complex<data_t> e2=exp(i*k*(r+s)*z);
    std::complex<data_t> e3=exp(i*k*(r-s)*z);
    std::complex<data_t> e4=exp(-i*k*(r-s)*z);
    std::complex<data_t> e5=exp(-i*k*(r+s)*z);
    M[0] =-2.0*r;   M[1] =1.0+rs;            M[2] =1.0-rs;              M[3] =M[2];         M[4] =M[1];     M[5] =-2.0*s;
    M[6] =4.0*mr;   M[7] =-mt-2.0*mrs;       M[8] =-mt+2.0*mrs;         M[9] =M[8];         M[10]=M[7];     M[11]=2.0*mt*s;
    M[12]=0;        M[13]=-2.0*ms+mt*s;      M[14]=2.0*ms-mt*s;         M[15]=M[13];        M[16]=M[14];    M[17]=0;
    M[18]=0;        M[19]=-mt*r+2.0*mr;      M[20]=M[19];               M[21]=mt*r-2.0*mr;  M[22]=M[21];    M[23]=0;
    M[24]=2.0*mt*r;   M[25]=-2.0*mrs-mt;       M[26]=2.0*mrs-mt;          M[27]=M[26];        M[28]=M[25];    M[29]=4.0*ms;
    M[30]=-4.0*mt*mr; M[31]=4.0*mr*ms+mt*mt;   M[32]=-4.0*mr*ms+mt*mt;    M[33]=M[32];        M[34]=M[31];    M[35]=-4.0*mt*ms;
    for (int j=0; j<6; j++){
        M[j*6+1] *= e2;
        M[j*6+2] *= e3;
        M[j*6+3] *= e4;
        M[j*6+4] *= e5;
    }
}

void compound_inv_field_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    data_t c=2*M_PI*f/k;
    if ((c == a) || (c == b)) c += EPS;
    std::complex<data_t> r=sqrt((std::complex<data_t>)((c/a)*(c/a)-1));
    std::complex<data_t> s=sqrt((std::complex<data_t>)((c/b)*(c/b)-1));
    data_t g=(b/c)*(b/c);
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    std::complex<data_t> e2=exp(-i*k*(r+s)*z);
    std::complex<data_t> e3=exp(-i*k*(r-s)*z);
    std::complex<data_t> e4=exp(i*k*(r-s)*z);
    std::complex<data_t> e5=exp(i*k*(r+s)*z);
    M[0] =g*(2*g-1)/r;                  M[1] =g*g/mr;                         M[2] =0;              M[3] =0;            M[4] =M[0]/(2*mu);  M[5] =M[1]/(2*mu);
    M[6] =g*g+(2*g-1)*(2*g-1)/(4.0*rs); M[7] =g*g/(2*mu)+g*(2*g-1)/(4.0*mrs); M[8] =-g/(4.0*ms);    M[9] =g/(4.0*mr);   M[10]=M[7];         M[11]=(g/(2*mu))*(g/(2*mu))*(1.0+1.0/rs);
    M[12]=2*g*g-M[6];                   M[13]=g*g/mu-M[7];                    M[14]=-M[8];          M[15]=M[9];         M[16]=M[13];        M[17]=(g/(2*mu))*(g/(2*mu))*(1.0-1.0/rs);
    M[18]=M[12];                        M[19]=M[13];                          M[20]=M[8];           M[21]=-M[9];        M[22]=M[13];        M[23]=M[17];
    M[24]=M[6];                         M[25]=M[7];                           M[26]=-M[8];          M[27]=-M[9];        M[28]=M[7];         M[29]=M[11];
    M[30]=g*(2*g-1)/s;                  M[31]=M[30]/(2*mu);                   M[32]=0;              M[33]=0;            M[34]=g*g/ms;       M[35]=M[34]/(2*mu);
    for (int j=0; j<6; j++){
        M[6+j] *= e2;
        M[12+j] *= e3;
        M[18+j] *= e4;
        M[24+j] *= e5;
    }
}

void compound_inv_field_matrix_c(std::complex<data_t> * M, const data_t f, const std::complex<data_t> k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    std::complex<data_t> c=2*M_PI*f/k;
    if ((c.real() == a) || (c.real() == b)) c += EPS;
    std::complex<data_t> r=sqrt_mb((c/a)*(c/a)-1.0);
    std::complex<data_t> s=sqrt_mb((c/b)*(c/b)-1.0);
    //r *= (data_t)sign(r.imag());
    //r *= (data_t)sign(k.real()*r.real()-k.imag()*r.imag());
    //s *= (data_t)sign(s.imag());
    //s *= (data_t)sign(k.real()*s.real()-k.imag()*s.imag());
    std::complex<data_t> g=(b/c)*(b/c);
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    std::complex<data_t> e2=exp(-i*k*(r+s)*z);
    std::complex<data_t> e3=exp(-i*k*(r-s)*z);
    std::complex<data_t> e4=exp(i*k*(r-s)*z);
    std::complex<data_t> e5=exp(i*k*(r+s)*z);
    M[0] =g*(2.0*g-1.0)/r;                  M[1] =g*g/mr;                         M[2] =0;              M[3] =0;            M[4] =M[0]/(2*mu);  M[5] =M[1]/(2*mu);
    M[6] =g*g+(2.0*g-1.0)*(2.0*g-1.0)/(4.0*rs); M[7] =g*g/(2*mu)+g*(2.0*g-1.0)/(4.0*mrs); M[8] =-g/(4.0*ms);    M[9] =g/(4.0*mr);   M[10]=M[7];         M[11]=(g/(2*mu))*(g/(2*mu))*(1.0+1.0/rs);
    M[12]=2.0*g*g-M[6];                   M[13]=g*g/mu-M[7];                    M[14]=-M[8];          M[15]=M[9];         M[16]=M[13];        M[17]=(g/(2*mu))*(g/(2*mu))*(1.0-1.0/rs);
    M[18]=M[12];                        M[19]=M[13];                          M[20]=M[8];           M[21]=-M[9];        M[22]=M[13];        M[23]=M[17];
    M[24]=M[6];                         M[25]=M[7];                           M[26]=-M[8];          M[27]=-M[9];        M[28]=M[7];         M[29]=M[11];
    M[30]=g*(2.0*g-1.0)/s;                  M[31]=M[30]/(2*mu);                   M[32]=0;              M[33]=0;            M[34]=g*g/ms;       M[35]=M[34]/(2*mu);
    for (int j=0; j<6; j++){
        M[6+j] *= e2;
        M[12+j] *= e3;
        M[18+j] *= e4;
        M[24+j] *= e5;
    }
}

void compound_propagator_matrix(std::complex<data_t> * M, const data_t f, const data_t k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    data_t c=2*M_PI*f/k;
    if ((c == a) || (c == b)) c += EPS;
    std::complex<data_t> r=sqrt((std::complex<data_t>)((c/a)*(c/a)-1));
    std::complex<data_t> s=sqrt((std::complex<data_t>)((c/b)*(c/b)-1));
    data_t g=(b/c)*(b/c);
    data_t t=2-(c/b)*(c/b);
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    std::complex<data_t> Ca=cos(k*r*z);
    std::complex<data_t> Cb=cos(k*s*z);
    std::complex<data_t> Sa=sin(k*r*z);
    std::complex<data_t> Sb=sin(k*s*z);
    std::complex<data_t> CaCb=Ca*Cb;
    std::complex<data_t> SaSb=Sa*Sb;
    std::complex<data_t> CaSb=Ca*Sb;
    std::complex<data_t> SaCb=Sa*Cb;
    M[0]=-g*g*(4*t-(t*t+4)*CaCb+(4.0*rs+t*t/rs)*SaSb);
    M[1]=-g*g/mu*((2+t)*(1.0-CaCb)+(2.0*rs+t/rs)*SaSb);
    M[2]=-i*g/mu*(CaSb/s+r*SaCb);
    M[3]=i*g/mu*(SaCb/r+s*CaSb);
    M[4]=M[1];
    M[5]=-(g/mu)*(g/mu)*(2.0+(1.0/rs+rs)*SaSb-2.0*CaCb);
    M[6]=g*g*mu*(2*t*(2+t)*(1.0-CaCb)+(8.0*rs+t*t*t/rs)*SaSb);
    M[7]=1.0+CaCb-M[0];
    M[8]=i*g*(t/s*CaSb+2.0*r*SaCb);
    M[9]=-i*g*(t/r*SaCb+2.0*s*CaSb);
    M[10]=M[7]-1.0;
    M[11]=-M[1];
    M[12]=-i*g*mu*(4.0*s*CaSb+t*t/r*SaCb);
    M[13]=M[9];
    M[14]=CaCb;
    M[15]=s/r*SaSb;
    M[16]=M[9];
    M[17]=-M[3];
    M[18]=i*g*mu*(4.0*r*SaCb+t*t/s*CaSb);
    M[19]=M[8];
    M[20]=r/s*SaSb;
    M[21]=M[14];
    M[22]=M[8];
    M[23]=-M[2];
    M[24]=M[6];
    M[25]=M[10];
    M[26]=M[8];
    M[27]=M[9];
    M[28]=M[7];
    M[29]=-M[1];
    M[30]=-(g*mu)*(g*mu)*(8*t*t*(1.0-CaCb)+(16.0*rs+t*t*t*t/rs)*SaSb);
    M[31]=-M[6];
    M[32]=-M[18];
    M[33]=-M[12];
    M[34]=-M[6];
    M[35]=M[0];
}

void compound_propagator_matrix_c(std::complex<data_t> * M, const data_t f, const std::complex<data_t> k, const data_t a, const data_t b, const data_t rho, const data_t z){
    std::complex<data_t> i(0,1);
    data_t mu=rho*b*b;
    std::complex<data_t> c=2*M_PI*f/k;
    if ((c.real() == a) || (c.real() == b)) c += EPS;
    std::complex<data_t> r=sqrt_mb((c/a)*(c/a)-1.0);
    std::complex<data_t> s=sqrt_mb((c/b)*(c/b)-1.0);
    // r *= (data_t)sign(r.imag());
    // r *= (data_t)sign(k.real()*r.real()-k.imag()*r.imag());
    // s *= (data_t)sign(s.imag());
    // s *= (data_t)sign(k.real()*s.real()-k.imag()*s.imag());
    std::complex<data_t> g=(b/c)*(b/c);
    std::complex<data_t> t=2.0-(c/b)*(c/b);
    std::complex<data_t> rs=r*s;
    std::complex<data_t> mrs=mu*r*s;
    std::complex<data_t> mr=mu*r;
    std::complex<data_t> ms=mu*s;
    std::complex<data_t> Ca=cos(k*r*z);
    std::complex<data_t> Cb=cos(k*s*z);
    std::complex<data_t> Sa=sin(k*r*z);
    std::complex<data_t> Sb=sin(k*s*z);
    std::complex<data_t> CaCb=Ca*Cb;
    std::complex<data_t> SaSb=Sa*Sb;
    std::complex<data_t> CaSb=Ca*Sb;
    std::complex<data_t> SaCb=Sa*Cb;
    M[0]=-g*g*(4.0*t-(t*t+4.0)*CaCb+(4.0*rs+t*t/rs)*SaSb);
    M[1]=-g*g/mu*((2.0+t)*(1.0-CaCb)+(2.0*rs+t/rs)*SaSb);
    M[2]=-i*g/mu*(CaSb/s+r*SaCb);
    M[3]=i*g/mu*(SaCb/r+s*CaSb);
    M[4]=M[1];
    M[5]=-(g/mu)*(g/mu)*(2.0+(1.0/rs+rs)*SaSb-2.0*CaCb);
    M[6]=g*g*mu*(2.0*t*(2.0+t)*(1.0-CaCb)+(8.0*rs+t*t*t/rs)*SaSb);
    M[7]=1.0+CaCb-M[0];
    M[8]=i*g*(t/s*CaSb+2.0*r*SaCb);
    M[9]=-i*g*(t/r*SaCb+2.0*s*CaSb);
    M[10]=M[7]-1.0;
    M[11]=-M[1];
    M[12]=-i*g*mu*(4.0*s*CaSb+t*t/r*SaCb);
    M[13]=M[9];
    M[14]=CaCb;
    M[15]=s/r*SaSb;
    M[16]=M[9];
    M[17]=-M[3];
    M[18]=i*g*mu*(4.0*r*SaCb+t*t/s*CaSb);
    M[19]=M[8];
    M[20]=r/s*SaSb;
    M[21]=M[14];
    M[22]=M[8];
    M[23]=-M[2];
    M[24]=M[6];
    M[25]=M[10];
    M[26]=M[8];
    M[27]=M[9];
    M[28]=M[7];
    M[29]=-M[1];
    M[30]=-(g*mu)*(g*mu)*(8.0*t*t*(1.0-CaCb)+(16.0*rs+t*t*t*t/rs)*SaSb);
    M[31]=-M[6];
    M[32]=-M[18];
    M[33]=-M[12];
    M[34]=-M[6];
    M[35]=M[0];
}

std::complex<data_t> free_solid_free(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){

    std::complex<data_t> answer;

    // propagator matrix
    std::vector<std::complex<data_t> > T (36);
    compound_propagator_matrix(T.data(),f, k, a[0], b[0], rho[0], d[0]);

    // One layer only
    if (n==1){
        answer = T[30];
        return answer;
    }

    else{
        std::vector<std::complex<data_t> > temp; temp = {T[0],T[6],T[12],T[18],T[24],T[30]};
        std::vector<std::complex<data_t> > v(6);
        for (int i=1; i<n; i++){

            // propagator matrix
            compound_propagator_matrix(T.data(), f, k, a[i], b[i], rho[i], d[i]);

            // v = T*v  v0 = T(:,1)
            #pragma omp parallel for
            for (int j=0; j<6; j++){
                v[j] = 0;
                for (int k=0; k<6; k++){
                    v[j] += T[j*6+k]*temp[k];
                }
            }
            temp = v;
        }

        answer = v[5];
        return answer;
    }

}

std::complex<data_t> rigid_solid_rigid(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    
    std::complex<data_t> answer;

    // propagator matrix
    std::vector<std::complex<data_t> > T (36);
    compound_propagator_matrix(T.data(),f, k, a[0], b[0], rho[0], d[0]);

    // One layer only
    if (n==1){
        answer = T[5];
        return answer;
    }

    else{
        std::vector<std::complex<data_t> > temp; temp = {T[5],T[11],T[17],T[23],T[29],T[35]};
        std::vector<std::complex<data_t> > v(6);
        for (int i=1; i<n; i++){

            // propagator matrix
            compound_propagator_matrix(T.data(), f, k, a[i], b[i], rho[i], d[i]);

            // v = T*v  v0 = T(:,6)
            #pragma omp parallel for
            for (int j=0; j<6; j++){
                v[j] = 0;
                for (int k=0; k<6; k++){
                    v[j] += T[j*6+k]*temp[k];
                }
            }
            temp = v;
        }

        answer = v[0];
        return answer;
    }

}

std::complex<data_t> free_solid_halfspace(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    
    // Field matrix (Q matrix) for the half space
    std::vector<std::complex<data_t> > Q (36);
    compound_field_matrix(Q.data(),f, k, a[n-1], b[n-1], rho[n-1], 0);

    std::complex<data_t> answer;
    // Only half space
    if (n==1){
        answer = Q[31];
        return answer;
    }

    // Single layer overlying a half space
    else if (n==2){
        // propagator matrix
        std::vector<std::complex<data_t> > T (36);
        compound_propagator_matrix(T.data(), f, k, a[0], b[0], rho[0], -d[0]);

        // T(6,:)*Q(:,2)
        std::complex<data_t> temp (0,0);
        for (int i=0; i<6; i++) temp += T[30+i]*Q[i*6+1];
        answer = temp;
        return answer;
    }

    // At least 2 layers overlying a half space 
    else {
        std::vector<std::complex<data_t> > T (36);
        std::vector<std::complex<data_t> > temp; temp = {Q[1],Q[7],Q[13],Q[19],Q[25],Q[31]};
        std::vector<std::complex<data_t> > v(6);
        for (int i=n-2; i>0; i--){
            // propagator matrix
            compound_propagator_matrix(T.data(), f, k, a[i], b[i], rho[i], -d[i]);

            // v = T*v  v0 = Q(:,2)
            #pragma omp parallel for
            for (int j=0; j<6; j++){
                v[j] = 0;
                for (int k=0; k<6; k++){
                    v[j] += T[j*6+k]*temp[k];
                }
            }
            temp = v;
        }
        // T(6,:)*v
        compound_propagator_matrix(T.data(), f, k, a[0], b[0], rho[0], -d[0]);
        std::complex<data_t> temp2 (0,0);
        for (int i=0; i<6; i++) temp2 += T[30+i]*v[i];
        answer = temp2;
        return answer;
    }
}

std::complex<data_t> halfspace_solid_halfspace(const data_t f, const data_t k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    
    if (n<2){
        fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
        throw std::logic_error("The number of layers must be at least 2.\n");
    }

    // Field matrix (Q matrix) for the top half space
    std::vector<std::complex<data_t> > Q (36);
    compound_field_matrix(Q.data(),f, k, a[0], b[0], rho[0], 0);

    std::complex<data_t> answer;
    // Adjacent two half spaces (Stoneley mode)
    if (n==2){
        // Inverse field matrix (invQ matrix) for the lower half space
        std::vector<std::complex<data_t> > invQ (36);
        compound_inv_field_matrix(invQ.data(),f, k, a[1], b[1], rho[1], 0);

        // invQ(5,:)*Q(:,5)
        std::complex<data_t> temp(0,0);
        for (int i=0; i<6; i++) temp += invQ[24+i]*Q[i*6+4];
        answer = temp;
        return answer;
    }

    // At least one layer embedded between the two halfspaces
    else {
        std::vector<std::complex<data_t> > T (36);
        std::vector<std::complex<data_t> > temp; temp = {Q[4],Q[10],Q[16],Q[22],Q[28],Q[34]};
        std::vector<std::complex<data_t> > v(6);
        data_t Z=0;
        for (int i=1; i<n-1; i++){
            // propagator matrix
            compound_propagator_matrix(T.data(), f, k, a[i], b[i], rho[i], d[i]);
            Z += d[i];

            // v = T*v  v0 = Q(:,5)
            #pragma omp parallel for
            for (int j=0; j<6; j++){
                v[j] = 0;
                for (int k=0; k<6; k++){
                    v[j] += T[j*6+k]*temp[k];
                }
            }
            temp = v;
        }
        // Inverse field matrix at depth Z for the lower halfspace
        compound_inv_field_matrix(T.data(), f, k, a[n-1], b[n-1], rho[n-1], Z);

        // invD(5,:)*v
        std::complex<data_t> temp2 (0,0);
        for (int i=0; i<6; i++) temp2 += T[24+i]*v[i];
        answer = temp2;
        return answer;
    }
}

std::complex<data_t> halfspace_solid_halfspace_leaky(const data_t f, const std::complex<data_t> k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    
    if (n<2){
        fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
        throw std::logic_error("The number of layers must be at least 2.\n");
    }

    // Field matrix (Q matrix) for the top half space
    std::vector<std::complex<data_t> > Q (36);
    // compound_field_matrix_c(Q.data(),f, k, a[0], b[0], rho[0], 0);
    compound_field_matrix(Q.data(),f, k.real(), a[0], b[0], rho[0], 0);

    std::complex<data_t> answer;
    // Adjacent two half spaces (Stoneley mode)
    if (n==2){
        // Inverse field matrix (invQ matrix) for the lower half space
        std::vector<std::complex<data_t> > invQ (36);
        // compound_inv_field_matrix_c(invQ.data(),f, k, a[1], b[1], rho[1], 0);
        compound_inv_field_matrix(invQ.data(),f, k.real(), a[1], b[1], rho[1], 0);

        // invQ(5,:)*Q(:,5)
        std::complex<data_t> temp(0,0);
        for (int i=0; i<6; i++) temp += invQ[24+i]*Q[i*6+4];
        answer = temp;
        return answer;
    }

    // At least one layer embedded between the two halfspaces
    else {
        std::vector<std::complex<data_t> > T (36);
        std::vector<std::complex<data_t> > temp; temp = {Q[4],Q[10],Q[16],Q[22],Q[28],Q[34]};
        std::vector<std::complex<data_t> > v(6);
        data_t Z=0;
        for (int i=1; i<n-1; i++){
            // propagator matrix
            compound_propagator_matrix_c(T.data(), f, k, a[i], b[i], rho[i], d[i]);
            // compound_propagator_matrix(T.data(), f, k.real(), a[i], b[i], rho[i], d[i]);
            Z += d[i];

            // v = T*v  v0 = Q(:,5)
            #pragma omp parallel for
            for (int j=0; j<6; j++){
                v[j] = 0;
                for (int k=0; k<6; k++){
                    v[j] += T[j*6+k]*temp[k];
                }
            }
            temp = v;
        }
        // Inverse field matrix at depth Z for the lower halfspace
        // compound_inv_field_matrix_c(T.data(), f, k, a[n-1], b[n-1], rho[n-1], Z);
        compound_inv_field_matrix(T.data(), f, k.real(), a[n-1], b[n-1], rho[n-1], Z);

        // invD(5,:)*v
        std::complex<data_t> temp2 (0,0);
        for (int i=0; i<6; i++) temp2 += T[24+i]*v[i];
        answer = temp2;
        return answer;
    }
}

std::complex<data_t> free_solid_halfspace_leaky(const data_t f, const std::complex<data_t> k, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    
    // Field matrix (Q matrix) for the half space
    std::vector<std::complex<data_t> > Q (36);
    // compound_field_matrix_c(Q.data(),f, k, a[n-1], b[n-1], rho[n-1], 0);
    compound_field_matrix(Q.data(),f, k.real(), a[n-1], b[n-1], rho[n-1], 0);

    std::complex<data_t> answer;
    // Only half space
    if (n==1){
        answer = Q[31];
        return answer;
    }

    // Single layer overlying a half space
    else if (n==2){
        // propagator matrix
        std::vector<std::complex<data_t> > T (36);
        compound_propagator_matrix_c(T.data(), f, k, a[0], b[0], rho[0], -d[0]);

        // T(6,:)*Q(:,2)
        std::complex<data_t> temp (0,0);
        for (int i=0; i<6; i++) temp += T[30+i]*Q[i*6+1];
        answer = temp;
        return answer;
    }

    // At least 2 layers overlying a half space 
    else {
        std::vector<std::complex<data_t> > T (36);
        std::vector<std::complex<data_t> > temp; temp = {Q[1],Q[7],Q[13],Q[19],Q[25],Q[31]};
        std::vector<std::complex<data_t> > v(6);
        for (int i=n-2; i>0; i--){
            // propagator matrix
            compound_propagator_matrix_c(T.data(), f, k, a[i], b[i], rho[i], -d[i]);

            // v = T*v  v0 = Q(:,2)
            #pragma omp parallel for
            for (int j=0; j<6; j++){
                v[j] = 0;
                for (int k=0; k<6; k++){
                    v[j] += T[j*6+k]*temp[k];
                }
            }
            temp = v;
        }
        // T(6,:)*v
        compound_propagator_matrix_c(T.data(), f, k, a[0], b[0], rho[0], -d[0]);
        std::complex<data_t> temp2 (0,0);
        for (int i=0; i<6; i++) temp2 += T[30+i]*v[i];
        answer = temp2;
        return answer;
    }
}

void free_solid_free_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    #pragma omp parallel for
    for (int i=0; i<nf; i++){
        for (int j=0; j<nk; j++){
            val[i*nk+j] = free_solid_free(f[i],k[j],n,a,b,rho,d);
        }
    }
}

void rigid_solid_rigid_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    #pragma omp parallel for
    for (int i=0; i<nf; i++){
        for (int j=0; j<nk; j++){
            val[i*nk+j] = rigid_solid_rigid(f[i],k[j],n,a,b,rho,d);
        }
    }
}

void free_solid_halfspace_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    #pragma omp parallel for
    for (int i=0; i<nf; i++){
        for (int j=0; j<nk; j++){
            val[i*nk+j] = free_solid_halfspace(f[i],k[j],n,a,b,rho,d);
        }
    }
}

void halfspace_solid_halfspace_map(std::complex<data_t> * val, const data_t* f, const int nf, const data_t* k, const int nk, const int n, const data_t* a, const data_t* b, const data_t* rho, const data_t* d){
    #pragma omp parallel for
    for (int i=0; i<nf; i++){
        for (int j=0; j<nk; j++){
            val[i*nk+j] = halfspace_solid_halfspace(f[i],k[j],n,a,b,rho,d);
        }
    }
}

std::complex<data_t> sqrt_mb(std::complex<data_t> x){
    std::complex<data_t> answer(0,0);
    data_t xr = x.real();
    data_t xi = x.imag();
    if (xi == 0){
        if (xr >= 0) answer.real(sqrt(xr));
        else answer.imag(-sqrt(-xr));
    }
    else {
        answer.real(sqrt((xr+sqrt(xr*xr+xi*xi))/2));
        answer.imag(xi/(2*answer.real()));
    }
    //answer = sqrt(x);
    return answer;
}

void print(std::complex<data_t> * M, int m, int n){
    for (int j=0; j<n; j++){
        for (int i=0; i<m; i++) fprintf(stderr,"%.5f + %.5fi\t",M[j*m+i].real(),M[j*m+i].imag());
        fprintf(stderr,"\n");
    }
}