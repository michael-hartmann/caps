
#ifndef CQUADPACK_EXPORT_H
#define CQUADPACK_EXPORT_H

#ifdef CQUADPACK_STATIC_DEFINE
#  define CQUADPACK_EXPORT
#  define CQUADPACK_NO_EXPORT
#else
#  ifndef CQUADPACK_EXPORT
#    ifdef cquadpack_EXPORTS
        /* We are building this library */
#      define CQUADPACK_EXPORT 
#    else
        /* We are using this library */
#      define CQUADPACK_EXPORT 
#    endif
#  endif

#  ifndef CQUADPACK_NO_EXPORT
#    define CQUADPACK_NO_EXPORT 
#  endif
#endif

#ifndef CQUADPACK_DEPRECATED
#  define CQUADPACK_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CQUADPACK_DEPRECATED_EXPORT
#  define CQUADPACK_DEPRECATED_EXPORT CQUADPACK_EXPORT CQUADPACK_DEPRECATED
#endif

#ifndef CQUADPACK_DEPRECATED_NO_EXPORT
#  define CQUADPACK_DEPRECATED_NO_EXPORT CQUADPACK_NO_EXPORT CQUADPACK_DEPRECATED
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define CQUADPACK_NO_DEPRECATED
#endif

#endif
