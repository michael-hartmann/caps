/**
 * @file   material.c
 * @author Michael Hartmann <caps@speicherleck.de>
 * @date   January, 2019
 * @brief  support for arbitrary dielectric functions
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "material.h"
#include "utils.h"
#include "misc.h"

/** @brief Helper function to parse strings
 *
 * Parse a string in the form of "key separator value" where key and value
 * represent floating numbers. If key or separator is not found, false is
 * returned. If the string is matched successfully, value is set.
 *
 * @param [in]  line string to parse
 * @param [in]  key key
 * @param [in]  separator separator
 * @param [out] value numerical value of the string "value"
 * @retval true parsing successful
 * @retval false parsing not successful
 */
static bool _parse(const char *line, const char *key, const char separator, double *value)
{
    char *p = strstr(line, key);
    if(p != NULL)
    {
        p = strchr(p+strlen(key), separator); /* find separator */
        if(p != NULL)
        {
            *value = strtodouble(p+1);
            return true;
        }
    }

    return false;
}

/** @brief Initialize material
 *
 * The material properties are read from the file given by filename.
 *
 * This function temporarily overwrites the value of LC_NUMERIC in the
 * environment. LC_NUMERIC is restored before returning from the function.
 *
 * Be aware that this function does not check every corner case, so it is
 * dangerous to read untrusted files.
 *
 * @param [in] filename path to material specification
 * @param [in] calL \f$L+R\f$, separation between plane and center of sphere
 * @retval material if successful
 * @retval NULL if file cannot be read or is in wrong format
 */
material_t *material_init(const char *filename, double calL)
{
    size_t points = 0;
    size_t size = 128; /* number of elems of array */
    char line[2048];
    FILE *f = fopen(filename, "r");

    if(f == NULL)
        return NULL;

    material_t *material = xmalloc(sizeof(material_t));

    material->calL = calL;

    memset(material->filename, '\0', sizeof(material->filename));
    strncpy(material->filename, filename, sizeof(material->filename)/sizeof(char)-sizeof(char));

    /* initialize to 0 */
    material->omegap_low  = 0;
    material->gamma_low   = 0;
    material->omegap_high = 0;
    material->gamma_high  = 0;

    material->xi    = xmalloc(size*sizeof(double));
    material->epsm1 = xmalloc(size*sizeof(double));

    while(fgets(line, sizeof(line)/sizeof(char), f) != NULL)
    {
        /* remove whitespace at beginning and end of line */
        strim(line);

        /* replace each tab by a single space */
        strrep(line, '\t', ' ');

        if(line[0] == '#')
        {
            _parse(line, "omegap_low",  '=', &material->omegap_low);
            _parse(line, "gamma_low",   '=', &material->gamma_low);
            _parse(line, "omegap_high", '=', &material->omegap_high);
            _parse(line, "gamma_high",  '=', &material->gamma_high);

            continue;
        }

        const char *p = strchr(line, ' ');
        if(p != NULL)
        {
            double xi = strtodouble(line), epsm1 = strtodouble(p)-1;

            /* increase buffer */
            if(points >= size)
            {
                size *= 2;
                material->xi    = xrealloc(material->xi,    size*sizeof(double));
                material->epsm1 = xrealloc(material->epsm1, size*sizeof(double));
            }

            material->xi[points]    = xi;
            material->epsm1[points] = epsm1;

            /* check if the values of xi are sorted (ascending) */
            if(points > 0 && xi <= material->xi[points-1])
            {
                material_free(material);
                material = NULL;
                goto out;
            }

            points++;
        }
    }

    material->xi_min = material->xi[0];
    material->xi_max = material->xi[points-1];
    material->points = points;

    /* convert from eV to rad/s */
    material->omegap_low  /= CAPS_HBAR_EV;
    material->omegap_high /= CAPS_HBAR_EV;
    material->gamma_low   /= CAPS_HBAR_EV;
    material->gamma_high  /= CAPS_HBAR_EV;

out:
    if(f != NULL)
        fclose(f);

    return material;
}

/** @brief Get extrapolation parameters
 *
 * For frequencies where there is no tabulated data available, the value of the
 * dielectric function will be extrapolated assuming Drude behaviour:
 * \f[
 *      \epsilon(\mathrm{i}\xi) = 1+\frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}
 * \f]
 *
 * The parameters for the plasma frequency \f$\omega_\mathrm{P}\f$ and the
 * relaxation frequency \f$\gamma\f$ for \f$\xi>\xi_\mathrm{max}\f$ and
 * \f$\xi<\xi_\mathrm{min}\f$ will be stored into omegap_high, gamma_high, and
 * omegap_low, gamma_low. If a pointer is NULL, the memory is not referenced.
 *
 * @param [in]  material material object
 * @param [out] omegap_low plasma frequency for high-frequency extrapolation (in rad/s)
 * @param [out] gamma_low relaxation frequency for high-frequency extrapolation (in rad/s)
 * @param [out] omegap_high plasma frequency for low-frequency extrapolation (in rad/s)
 * @param [out] gamma_high relaxation frequency for low-frequency extrapolation (in rad/s)
 */
void material_get_extrapolation(material_t *material, double *omegap_low, double *gamma_low, double *omegap_high, double *gamma_high)
{
    if(omegap_low != NULL)
        *omegap_low = material->omegap_low;
    if(gamma_low != NULL)
        *gamma_low = material->gamma_low;
    if(omegap_high != NULL)
        *omegap_high = material->omegap_high;
    if(gamma_high != NULL)
        *gamma_high = material->gamma_high;
}

/** @brief Free material object
 *
 * @param material material object
 */
void material_free(material_t *material)
{
    if(material != NULL)
    {
        xfree(material->xi);
        xfree(material->epsm1);
        xfree(material);
    }
}

/** @brief Print information about object to stream
 *
 * Print information (filename, number of points, \f$\xi_\mathrm{min}\f$,
 * \f$\xi_\mathrm{max}\f$, ...) to stream. If prefix is not NULL, each line
 * will start with the string given in prefix.
 *
 * @param [in] material material object
 * @param [in] stream output stream (e.g. stdout)
 * @param [in] prefix prefix for each line or NULL
 */
void material_info(material_t *material, FILE *stream, const char *prefix)
{
    if(prefix == NULL)
        prefix = "";

    fprintf(stream, "%sfilename    = %s\n",    prefix, material->filename);
    fprintf(stream, "%spoints      = %zu\n",   prefix, material->points);
    fprintf(stream, "%scalL        = %gm\n",   prefix, material->calL);
    fprintf(stream, "%sxi_min      = %g/m\n",  prefix, material->xi_min);
    fprintf(stream, "%sxi_max      = %g/m\n",  prefix, material->xi_max);
    fprintf(stream, "%somegap_high = %geV\n",  prefix, CAPS_HBAR_EV*material->omegap_high);
    fprintf(stream, "%sgamma_high  = %geV\n",  prefix, CAPS_HBAR_EV*material->gamma_high);
    fprintf(stream, "%sgamma_low   = %geV\n",  prefix, CAPS_HBAR_EV*material->gamma_low);
    fprintf(stream, "%somegap_low  = %geV\n",  prefix, CAPS_HBAR_EV*material->omegap_low);
}

/** @brief Dielectric function for material
 *
 * Return the dielectric function \f$\epsilon(\mathrm{i}\xi)-1\f$ for the
 * material. For frequencies greater (smaller) than the maximum (minimum)
 * tabulated frequency, an extrapolation using a Drude model is used. For the
 * tabulated values linear interpolation is used.
 *
 * @param [in] xi frequency in rad/s
 * @param [in] args material (must be of type material_t *)
 */
double material_epsilonm1(double xi, void *args)
{
    material_t *self = (material_t *)args;

    const double xi_min = self->xi_min, xi_max = self->xi_max;

    if(xi < xi_min)
    {
        /* in scaled units */
        double omegap = self->omegap_low;
        double gamma_ = self->gamma_low;

        return POW_2(omegap)/(xi*(xi+gamma_));
    }
    else if(xi > xi_max)
    {
        /* in scaled units */
        double omegap = self->omegap_high;
        double gamma_ = self->gamma_high;

        return POW_2(omegap)/(xi*(xi+gamma_));
    }

    /* do binary search */
    int left = 0, right = self->points-1;

    /* find indices of right and left value */
    while(right-left != 1)
    {
        int middle = (right+left)/2;
        double xi_middle = self->xi[middle];
        if(xi_middle > xi)
            right = middle;
        else
            left = middle;
    }

    const double xi_lower = self->xi[left], xi_upper = self->xi[right];
    const double epsm1_lower = self->epsm1[left], epsm1_upper = self->epsm1[right];

    /* linear interpolation */
    return epsm1_lower + (xi-xi_lower)*(epsm1_upper-epsm1_lower)/(xi_upper-xi_lower);
}
