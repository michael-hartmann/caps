/**
 * @file   utils.c
 * @author Michael Hartmann <caps@speicherleck.de>
 * @date   January, 2018
 * @brief  wrappers for malloc, calloc realloc, and a few more useful functions
 */

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#endif

#include "utils.h"

/** @brief Wrapper for malloc
 *
 * This function is a wrapper for malloc. If malloc fails \ref TERMINATE is
 * called.
 *
 * @param size size of bytes to allocate
 * @retval ptr pointer to memory
 */
void *xmalloc(size_t size)
{
    void *p = malloc(size);

    TERMINATE(p == NULL, "malloc failed: size=%zu", size);

    return p;
}


/** @brief Wrapper for calloc
 *
 * This function is a wrapper for calloc. If calloc fails \ref TERMINATE is
 * called.
 *
 * @param nmemb number of elements
 * @param size size of each element
 * @retval ptr pointer to memory
 */
void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);

    TERMINATE(p == NULL, "calloc failed: nmemb=%zu, size=%zu", nmemb, size);

    return p;
}


/** @brief Wrapper for realloc
 *
 * This function is a wrapper for realloc. If realloc fails \ref TERMINATE is
 * called.
 *
 * @param p ptr to old memory
 * @param size size
 * @retval newptr pointer to new memory
 */
void *xrealloc(void *p, size_t size)
{
    void *p2 = realloc(p, size);

    TERMINATE(p2 == NULL, "realloc failed: p=%p, size=%zu", p, size);

    return p2;
}

/** @brief Seconds since 01/01/1970
 *
 * This function returns the seconds since 1st Jan 1970 in µs precision.
 *
 * @retval time seconds since 1st Jan 1970
 */
double now(void)
{
#ifdef _WIN32
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME system_time;
    FILETIME file_time;
    uint64_t time;

    GetSystemTime(&system_time);
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    uint64_t sec  = (uint64_t)((time - EPOCH) / 10000000L);
    uint64_t usec = (uint64_t)(system_time.wMilliseconds * 1000);
    return sec + usec*1e-6;
}
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
#endif
}

/** @brief Disable buffering to stderr and stdout
 *
 */
void disable_buffering(void)
{
    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
}

/** @brief Replace character by different character in string
 *
 * Replace occurence of a by b in the string s.
 *
 * @param [in,out] s string, terminated by \0
 * @param [in] a    character to replace
 * @param [in] b    substitute
 */
void strrep(char *s, const char a, const char b)
{
    while(*s != '\0')
    {
        if(*s == a)
            *s = b;
        s++;
    }
}

/** @brief Check if strings are equal
 *
 * Check if the the strings s1 and s2 are equal (case insensitive).
 *
 * @param [in] s1   first string
 * @param [in] s2   second string
 */
int strcaseequal(const char *s1, const char *s2)
{
    if(strlen(s1) != strlen(s2))
        return 0;

    for(size_t i = 0; i < strlen(s1); i++)
        if(tolower(s1[i]) != tolower(s2[i]))
            return 0;

    return 1;
}

/** @brief Remove whitespace at beginng and end of string
 *
 * If str is NULL the function doesn't do anything. Otherwise, trailing
 * whitespace and whitespace at the beginning of the string are removed.
 *
 * @param str string
 */
void strim(char *str)
{
    if(str == NULL)
        return;

    /* trim left */
    size_t i, len = strlen(str);
    for(i = 0; i < len; i++)
        if(!isspace(str[i]))
            break;

    /* i is the number of whitespace characters at the beginning of str */
    if(i)
    {
        memmove(str, str+i, len-i);
        str[len-i] = '\0';
    }

    /* trim right */
    len = strlen(str);
    for(i = 0; i < len; i++)
        if(!isspace(str[len-i-1]))
            break;

    /* i is the number of whitespace characters at the end of str */
    str[len-i] = '\0';
}

/**
 * @brief Convert string to float
 *
 * Convert the string str to a double.
 *
 * @param [in]  str         string to parse
 * @retval      value       converted value
 */
double strtodouble(const char *str)
{
    double intpart = 0, fracpart = 0, exponent = 0;
    int sign = +1, len = 0, conversion = 0;

    // skip whitespace
    while(isspace(*str))
        str++;

    // check for sign (optional; either + or -)
    if(*str == '-')
    {
        sign = -1;
        str++;
    }
    else if(*str == '+')
        str++;

    // check for nan and inf
    if(tolower(str[0]) == 'n' && tolower(str[1]) == 'a' && tolower(str[2]) == 'n')
        return NAN;
    if(tolower(str[0]) == 'i' && tolower(str[1]) == 'n' && tolower(str[2]) == 'f')
        return sign*INFINITY;

    // find number of digits before decimal point
    {
        const char *p = str;
        len = 0;
        while(isdigit(*p))
        {
            p++;
            len++;
        }
    }

    if(len)
        conversion = 1;

    // convert intpart part of decimal point to a float
    {
        double f = 1;
        for(int i = 0; i < len; i++)
        {
            int v = str[len-1-i] - '0';
            intpart += v*f;
            f *= 10;
        }
        str += len;
    }

    // check for decimal point (optional)
    if(*str == '.')
    {
        const char *p = ++str;

        // find number of digits after decimal point
        len = 0;
        while(isdigit(*p))
        {
            p++;
            len++;
        }

        if(len)
            conversion = 1;

        // convert fracpart part of decimal point to a float
        double f = 0.1;
        for(int i = 0; i < len; i++)
        {
            int v = str[i] - '0';
            fracpart += v*f;
            f *= 0.1;
        }

        str = p;
    }

    if(conversion && (*str == 'e' || *str == 'E'))
    {
        int expsign = +1;
        const char *p = ++str;

        if(*p == '+')
            p++;
        else if(*p == '-')
        {
            expsign = -1;
            p++;
        }

        str = p;
        len = 0;
        while(isdigit(*p))
        {
            len++;
            p++;
        }

        int f = 1;
        for(int i = 0; i < len; i++)
        {
            int v = str[len-1-i]-'0';
            exponent += v*f;
            f *= 10;
        }

        exponent *= expsign;
    }

    if(!conversion)
        return NAN;

    return sign*(intpart+fracpart)*pow(10, exponent);
}
