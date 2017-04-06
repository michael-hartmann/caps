#ifndef MATERIAL_H
#define MATERIAL_H

typedef struct {
    char filename[512];
    double calL;
    double omegap_low, gamma_low;
    double omegap_high, gamma_high;
    double xi_min, xi_max;
    size_t points;
    double *xi, *epsm1;
} material_t;

material_t *material_init(const char *filename, double calL);
void material_info(material_t *material, FILE *stream, const char *prefix);
void material_free(material_t *material);

void material_get_extrapolation(material_t *material, double *omegap_low, double *gamma_low, double *omegap_high, double *gamma_high);

double material_epsilonm1(double xi_scaled, void *args);

#endif
