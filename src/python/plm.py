from matplotlib import pyplot as plt
from libcasimir import sfunc

l = sfunc.Plm(4,2,4.1,1)
for i,elem in enumerate(l):
    print(i,elem)

plt.plot(l)
plt.show()
