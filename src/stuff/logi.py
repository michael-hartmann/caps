from math import log

maximum = 20000

print("#include <math.h>")
print()
print("#include \"logi.h\"")
print()
print("static double lookup[] = { NAN, ")

for i in range(1,maximum):
    print("    %.19g, /* log(%d)*/" % (log(i),i))

print("    %.19g /* log(%d)*/" % (log(maximum),maximum))
print("};");
print()
print("static int elems = sizeof(lookup)/sizeof(lookup[0]);")
print()
print("double logi(int x)")
print("{")
print("    if(x < elems)")
print("        return lookup[x];")
print("    return log(x);")
print("}")
