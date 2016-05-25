import matplotlib.pylab as plab
import h5py
import numpy as np

def main():

    id_diss = ["/home/jorgehog/tmp/extraneighbor.h5", 1464174780766]
    id_growth = ["/tmp/extraneighbor.h5", 1464181228224]

    eq_times = []
    for id_set, name in zip([id_growth, id_diss], ["growth", "diss"]):
        path, id = id_set

        f = h5py.File(path, 'r')

        group = f["30x30"][str(id)]
        eq_set = group["eq_storedEventValues"][0]
        neq_set = group["omega_storedEventValues"][0]

        eq_time = group["eq_storedEventValues"][1]
        neq_time = group["omega_storedEventValues"][1]

        eq_times.append(eq_time[-1])
        full_set = np.concatenate([eq_set, neq_set], axis=0)
        full_time = np.concatenate([eq_time, neq_time + eq_time[-1] + neq_time[1]])

        np.save("/tmp/noneq_neigz_%s_cov.npy" % name, full_set)
        np.save("/tmp/noneq_neigz_%s_time.npy" % name, full_time)

    np.save("/tmp/noneq_neigz_eqtimes.npy", eq_times)

if __name__ == "__main__":
    main()