#ifndef MATERIAL_H
#define MATERIAL_H

#ifdef __cplusplus
extern "C" {
#endif

/** material_t data type */
typedef struct {
    char filename[512]; /**< material filename or \0\0... */
    double calL;        /**< \f$L+R\f$ */
    double xi_min;      /**< lower border of tabulated frequencies */
    double xi_max;      /**< upper border of tabulated frequencies */
    size_t points;      /**< number of points */
    double *xi;         /**< tabulated frequencies \f$\xi\f$ */
    double *epsm1;      /**< tabulated dielectric function, \f$\epsilon(\mathrm{i}\xi)-1\f$ */
    double omegap_low;  /**< plasma frequency for low frequency extrapolation */
    double gamma_low;   /**< relaxation frequency for low frequency extrapolation */
    double omegap_high; /**< plasma frequency for high frequency extrapolation */
    double gamma_high;  /**< relaxation frequency for hight frequency extrapolation */
} material_t;

material_t *material_init(const char *filename, double calL);
void material_info(material_t *material, FILE *stream, const char *prefix);
void material_free(material_t *material);

void material_get_extrapolation(material_t *material, double *omegap_low, double *gamma_low, double *omegap_high, double *gamma_high);

double material_epsilonm1(double xi, void *args);

#ifdef __cplusplus
}
#endif

#endif
