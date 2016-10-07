#ifndef __GAUNT_H
#define __GAUNT_H

inline double gaunt_log_a0(int n, int nu, int m);
inline double gaunt_a0(int n,int nu,int m);
void gaunt(const int n, const int nu, const int m, double a_tilde[]);

inline int gaunt_qmax(const int n, const int nu, const int m);

#endif
