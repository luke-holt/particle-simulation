#ifndef UTIL_H
#define UTIL_H

#define UTIL_ASSERT(c) \
    do { \
        if (!(c)) { \
            util_log("ASSERT", "%s:%d:%s assertion failed '%s'", \
                     __FILE__, __LINE__, __func__, #c); \
            abort(); \
        } \
    } while (0)

#define da_init(da) \
    do { \
        (da)->count = 0; (da)->capacity = 256; \
        (da)->items = (typeof((da)->items))malloc((da)->capacity*sizeof(*(da)->items)); \
        UTIL_ASSERT((da)->items); \
    } while (0)
#define da_delete(da) \
    do { \
        free((da)->items); (da)->items = NULL; \
        (da)->capacity = (da)->count = 0; \
    } while (0)
#define da_resize(da, cap) \
    do { \
        (da)->items = (typeof((da)->items))realloc((da)->items, sizeof(*(da)->items)*cap); \
        UTIL_ASSERT((da)->items); (da)->capacity = cap; \
    } while (0)
#define da_append(da, item) \
    do { \
        if ((da)->count >= (da)->capacity) \
            da_resize(da, (da)->capacity*2); \
        (da)->items[(da)->count++] = (item); \
    } while (0)


void util_log(const char *tag, const char *fmt, ...);

#endif // UTIL_H


// #define UTIL_IMPLEMENTATION
#ifdef UTIL_IMPLEMENTATION

#include <stdio.h>
#include <stdarg.h>

void
util_log(const char *tag, const char *fmt, ...)
{
    if (tag) fprintf(stdout, "[%s] ", tag);
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
}

#endif // UTIL_IMPLEMENTATION
