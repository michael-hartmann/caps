from glob import glob

for filename in glob("orig/*.txt"):
    pos1 = filename.find("/")
    pos2 = filename.find("nm")
    R = filename[pos1+1:pos2]

    with open("data_%snm.csv" % R, "w") as fout:
        print("# Casimir force F for a sphere of radius R above a plate at various separations L", file=fout)
        print("# measured by Ricardo Decca", file=fout)
        print("#", file=fout)
        print("# R = %snm" % R, file=fout)
        print("#", file=fout)
        print("# F_PR = -pi^3*hbar*c*R/(360*L^3)", file=fout)
        print("#", file=fout)
        print("# L (nm), F/F_PR", file=fout)
        with open(filename) as fin:
            for line in fin:
                line = line.strip()
                if len(line) == 0:
                    continue
                L, F = line.split()

                print("%s, %s" % (L, F), file=fout)
