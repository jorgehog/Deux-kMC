__author__ = 'jorgehog'

import os
import re
import h5py
import glob


class ParseKMCHDF5:

    def __init__(self, which_file):
        self.path, self.filename = os.path.split(which_file)

        self.files = []
        self.filenames = []

        match = re.findall("(.*)\_\d+\.h5", self.filename)

        if match:
            name = match[0]

            filenames = glob.glob1(self.path, name + "_*.h5")

            for file in filenames:

                fullpath = os.path.join(self.path, file)
                self.files.append(h5py.File(fullpath, 'r'))
                self.filenames.append(file)

        else:
            self.files = [h5py.File(which_file, 'r')]

        self.only_n = None

    def __iter__(self):

        for file in self.files:
            for l, run in file.items():
                L, W = [int(x) for x in re.findall("(\d+)x(\d+)", l)[0]]
                for run_id, data in run.items():

                    yield data, L, W, run_id


    def set_only_n(self, n):
        self.only_n = n

    def getfiles(self):
        return self.files

    def close(self):

        for file in self.files:
            file.close()

def getIgnisData(data, name):

    if not data.attrs["storeIgnisData"]:
        raise RuntimeError("No ignis data stored.")

    idx = list(data["ignisEventDescriptions"]).index(name)

    return data["ignisData"][idx]

if __name__ == "__main__":
    obj = ParseKMCHDF5("/home/jorgehog/code/build-kMC-Desktop_Qt_5_3_GCC_64bit-Release/apps/Quasi2DLoaded/Quasi2D.h5")

    for potential, alpha, mu, em, c, ignis_index_map, data, n in obj:
        print alpha, data["ignisData"][ignis_index_map["SurfaceSize"], -1]