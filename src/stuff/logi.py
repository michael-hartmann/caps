from math import log,lgamma

maximum = 20000

print("""
#include <math.h>

static double lookup_logi[] = {
    NAN,""")

for i in range(1,maximum):
    print("    %.20g, /* log(%d)*/" % (log(i),i))

print("    %.20g /* log(%d)*/" % (log(maximum),maximum))
print("};");

print("static double lookup_lfaci[] = {")
for i in range(0,maximum):
    print("    %.20g, /* lgamma(1+%d)*/" % (lgamma(1+i),i))

print("    %.20g /* lgamma(1+%d)*/" % (lgamma(1+maximum),maximum))
print("};");

print("""
static int elems_lookup_logi  = sizeof(lookup_logi) /sizeof(lookup_logi [0]);
static int elems_lookup_lfaci = sizeof(lookup_lfaci)/sizeof(lookup_lfaci[0]);

double logi(int x)
{
    if(x < elems_lookup_logi)
        return lookup_logi[x];
    return log(x);
}

double lfaci(int x)
{
    if(x < elems_lookup_lfaci)
        return lookup_lfaci[x];
    else
        return lgamma(1+x);
}

""")
