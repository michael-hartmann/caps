from libcasimir import Casimir
from mpi4py import MPI
from time import sleep


class CasimirMPI:
    def __init__(self, LbyR, ldim=None, tolerance=None, threshold=None):
        self.LbyR = LbyR
        self.ldim = ldim
        self.tolerance = tolerance
        self.threshold = None

        # mpi
        self.comm = MPI.COMM_WORLD
        self.size = comm.Get_size()
        self.idle = [True] * self.size

    def is_idle(self, slave):
        return self.idle[slave]

    def submit(self, mesg):
        comm = self.comm
        for slave in range(1,self.size):
            if self.is_idle(slave):
                comm.send(mesg, dest=slave)
                self.idle[slave] = False
                return True
        return False

    def retrieve(self):
        results = []
        comm = self.comm
        for slave in range(1,self.size):
            if self.idle[slave] == False:
                if comm.Iprobe(source=slave):
                    results.append(comm.recv(source=slave))
                    self.idle[slave] = True
        return results

    def get_running(self):
        return self.size-sum(self.idle[1:])

    def slave(self):
        comm = self.comm
        LbyR = self.LbyR
        ldim = self.ldim
        tolerance = self.tolerance
        threshold = self.threshold

        casimir = Casimir(LbyR, ldim=ldim, tolerance=tolerance, threshold=threshold)

        while True:
            request = comm.recv(source=0)
            if "exit" in d:
                exit(0)
            
            xi,m = request["xi"], request["m"]
            logdetD = casimir.logdetD(xi, m)

            answer = { "xi": xi, "m": m, "logdetD": logdetD }
            comm.send(answer, dest=0)


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        # master
        c = CasimirMPI()

        i = 0
        while i < 20:
            while c.submit(i):
                print("submitted", i)
                i=i+1

            print(c.retrieve(), c.idle)
            sleep(1)
        print("exit")
        comm.Abort()
        exit(0)
    else:
        slave()
