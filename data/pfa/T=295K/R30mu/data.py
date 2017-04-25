from glob import glob
from sys import path
path.append("..")
from casimir import print_force, print_gradient

def stepsize_force(L):
    if L > 300e-9:
        return 2
    return 1

def stepsize_gradient(L):
    if L > 500e-9:
        return 3
    if L > 300e-9:
        return 2
    return 1

if __name__ == "__main__":
    files = glob("gold_eta9/slurm-*.out")

    f_force    = open("force_30000.csv", "w")
    print_force(files, stepsize_force, f=f_force)

    f_gradient = open("gradient_30000.csv", "w")
    print_gradient(files, stepsize_gradient, f=f_gradient)
