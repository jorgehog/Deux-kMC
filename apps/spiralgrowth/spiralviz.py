#DCVIZ

from DCViz_sup import DCVizPlotter


class IgnisOut(DCVizPlotter):

    nametag = "ignisSOS"

    def plot(self, data):

        names = self.loader.get_metadata()[0].split()

        for column, name in zip(data, names):

            if name.startswith("AverageHeight"):
                self.subfigure.plot(column)
                break


