__author__ = 'jorgehog'

import os
import re
import h5py
import glob


class ParseKMCHDF5:

    def __init__(self, which_file):
        self.path, filename = os.path.split(which_file)

        self.files = []
        self.filenames = []

        match = re.findall("(.*)\_\d+\.h5", filename)

        if match:
            name = match[0]

            filenames = glob.glob1(self.path, name + "_*.h5")

            for file in filenames:

                fullpath = os.path.join(self.path, file)
                self.files.append(fullpath)
                self.filenames.append(file)

        else:
            self.files = [which_file]

        self.only_n = None

    def __iter__(self):

        for self.filename in self.files:

            file = h5py.File(self.filename, 'r')

            for l, run in file.items():
                L, W = [int(x) for x in re.findall("(\d+)x(\d+)", l)[0]]
                for run_id, data in run.items():
                    yield data, L, W, run_id

            file.close()

    @staticmethod
    def get_ignis_data(data, name):
        return data["ignisData"][[desc.split("@")[0] for desc in list(data["ignisEventDescriptions"])[0]].index(name)]

    def set_only_n(self, n):
        self.only_n = n

    def getfiles(self):
        return self.files


def getIgnisData(data, name):

    if not data.attrs["storeIgnisData"]:
        raise RuntimeError("No ignis data stored.")

    idx = list(data["ignisEventDescriptions"]).index(name)

    return data["ignisData"][idx]

if __name__ == "__main__":
    obj = ParseKMCHDF5("/home/jorgehog/code/build-kMC-Desktop_Qt_5_3_GCC_64bit-Release/apps/Quasi2DLoaded/Quasi2D.h5")

    for potential, alpha, mu, em, c, ignis_index_map, data, n in obj:
        print alpha, data["ignisData"][ignis_index_map["SurfaceSize"], -1]