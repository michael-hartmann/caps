from glob import glob
from sys import path
path.append("..")
from casimir import print_force, print_gradient

def stepsize_force(L):
    if L > 400e-9:
        return 3
    if L > 200e-9:
        return 2
    return 1

def stepsize_gradient(L):
    if L > 600e-9:
        return 5
    if L > 500e-9:
        return 4
    if L > 400e-9:
        return 3
    if L > 200e-9:
        return 2
    return 1

if __name__ == "__main__":
    files = glob("gold_eta8/slurm-*.out")

    f_force    = open("force_100000.csv", "w")
    print_force(files, stepsize_force, f=f_force)

    f_gradient = open("gradient_100000.csv", "w")
    print_gradient(files, stepsize_gradient, f=f_gradient)
