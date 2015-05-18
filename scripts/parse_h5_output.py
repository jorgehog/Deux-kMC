__author__ = 'jorgehog'

import os
import re
import h5py
import glob


class ParseKMCHDF5:

    def __init__(self, which_file):
        self.path, self.filename = os.path.split(which_file)

        self.file = None
        self.file_number = None
        self.files = []
        self.filenames = []

        self.skipWhen = lambda alpha, mu, E0, s0, r0, n: False

        match = re.findall("(.*)\_\d+\.h5", self.filename)

        if match:
            name = match[0]

            filenames = glob.glob1(self.path, name + "_*.h5")

            for file in filenames:

                fullpath = os.path.join(self.path, file)
                self.files.append(fullpath)
                self.filenames.append(file)

        else:
            self.files = [h5py.File(which_file, 'r')]

        self.only_n = None

    def __iter__(self):

        for self.file_number, self.filepath in enumerate(self.files):
            self.file = h5py.File(self.filepath, 'r')

            for l, run in self.file.items():
                L, W = [int(x) for x in re.findall("(\d+)x(\d+)", l)[0]]
                for potential in run.keys():

                    try:
                        if "_n_" in potential:
                            alpha, mu, E0, s0, r0, n = [float(re.findall("%s\_(-?\d+\.?\d*[eE]?-?\d*|-?nan)" % ID, potential)[0]) for ID in
                                                     ["alpha", "mu", "E0", "s0", "r0", "n"]]
                        else:
                            alpha, mu, E0, s0, r0 = [float(re.findall("%s\_(-?\d+\.?\d*[eE]?-?\d*|-?nan)" % ID, potential)[0]) for ID in
                                                     ["alpha", "mu", "E0", "s0", "r0"]]
                            n = 0
                    except:
                        raise ValueError("invalid potential: %s" % potential)

                    if self.skipWhen(alpha, mu, E0, s0, r0, n):
                        continue

                    data = run[potential]

                    if "nNeighbors" in data.attrs.keys():
                        neighbors = data.attrs["nNeighbors"]
                    else:
                        neighbors = 0

                    _ignis_index_map = {}

                    if data.attrs["storeIgnisData"]:
                        for i, name in enumerate(data["ignisEventDescriptions"][0]):
                            name_strip = str(name).split("@")[0]
                            _ignis_index_map[name_strip] = i

                    yield L, W, potential, alpha, mu, E0, s0, r0, neighbors, _ignis_index_map, data, n

            self.file.close()


        self.file = None
        self.file_number = None

    def get_data(self, file_number, name):

        if self.file_number != file_number:
            if self.file is not None:
                self.file.close()

            self.file = h5py.File(self.files[file_number], 'r')

        self.file_number = file_number
        data = self.file[name]

        return data

    def set_only_n(self, n):
        self.only_n = n

    def getfiles(self):
        return self.files

    def close(self):

        for file in self.files:
            file.close()

if __name__ == "__main__":
    obj = ParseKMCHDF5("/home/jorgehog/code/build-kMC-Desktop_Qt_5_3_GCC_64bit-Release/apps/Quasi2DLoaded/Quasi2D.h5")

    for potential, alpha, mu, em, c, ignis_index_map, data, n in obj:
        print alpha, data["ignisData"][ignis_index_map["SurfaceSize"], -1]