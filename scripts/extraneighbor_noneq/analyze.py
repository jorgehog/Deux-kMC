import h5py
import numpy as np

def main():
    id_diss_low = ["/home/jorgehog/Dropbox/PhD/papers/paper3/data/extran_noneq/extraneighbor_0.h5", 1464212767497]
    id_diss_high = ["/home/jorgehog/Dropbox/PhD/papers/paper3/data/extran_noneq/extraneighbor_0.h5", 1464217682529]
    id_growth_low = ["/home/jorgehog/Dropbox/PhD/papers/paper3/data/extran_noneq/extraneighbor_2.h5", 21464212760854]
    id_growth_high = ["/home/jorgehog/Dropbox/PhD/papers/paper3/data/extran_noneq/extraneighbor_1.h5", 11464217683030]

    ids = [[id_diss_high, id_growth_high],
           [id_diss_low, id_growth_low]]

    eq_times = [[], []]
    F0s = []

    L = 30
    W = 30
    g = "%dx%d" % (L, W)

    for i, id_set in enumerate(ids):

        for j, id_pair in enumerate(id_set):

            path, id = id_pair

            f = h5py.File(path, 'r')

            group = f[g][str(id)]

            if j == 0:
                F0s.append(group.attrs["F0"])

            eq_set = group["eq_storedEventValues"][0]
            neq_set = group["omega_storedEventValues"][0]

            eq_time = group["eq_storedEventValues"][1]
            neq_time = group["omega_storedEventValues"][1]

            eq_times[i].append(eq_time[-1])
            full_set = np.concatenate([eq_set, neq_set])
            full_time = np.concatenate([eq_time, neq_time + eq_time[-1] + neq_time[1]])

            np.save("/tmp/noneq_neigz_%d%d_cov.npy" % (i, j), full_set/float(L*W))
            np.save("/tmp/noneq_neigz_%d%d_time.npy" % (i, j), full_time)

    np.save("/tmp/noneq_neigz_eqtimes.npy", eq_times)
    np.save("/tmp/noneq_neigz_F0s.npy", F0s)

if __name__ == "__main__":
    main()