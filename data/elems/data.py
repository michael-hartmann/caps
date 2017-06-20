from os import system

cmd = "../../src/matrix_elems %g %g %d %d > %s"

eta = 5
m    = (1,10,)
nT   = (1, 5)
LbyR = (0.02,)

no = 0
for LbyR_ in LbyR:
    for m_ in m:
        for nT_ in nT:
            lmax = int(eta/LbyR_)

            filename = "%04d.csv" % no
            with open(filename, "w") as f:
                print("# LbyR=%g, m=%d, nT=%g, lmax=%d" % (LbyR_, m_, nT_, lmax), file=f)

            system(cmd % (LbyR_, nT_, m_, lmax, filename))

            no += 1
