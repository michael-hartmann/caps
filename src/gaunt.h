#ifndef __GAUNT_H
#define __GAUNT_H

double gaunt_log_a0(int n, int nu, int m);
double gaunt_a0(int n,int nu,int m);
void gaunt(const int n, const int nu, const int m, double a_tilde[]);

int gaunt_qmax(const int n, const int nu, const int m);

#endif
