import numpy as np

from pfs.ga.common.util import safe_deep_copy

class Associations():
    """
    Represents a list of fiber-target associations.
    """

    def __init__(self, assoc=None, orig=None):
        if not isinstance(orig, Associations):
            self.__assoc = assoc or []      # List of index arrays
        else:
            self.__assoc = assoc or safe_deep_copy(orig.__assoc)

    def __get_assoc(self):
        return self.__assoc

    assoc = property(__get_assoc)

    def find_multiplets(self):
        """
        Find targets which are associated with more than one fiber.
        """

        all = np.concatenate(self.__assoc)
        counts = np.bincount(all)
        ix = np.where(counts > 1)[0]
        return ix, counts[ix]


    def sort(self, sorting):
        """
        Sorts the associations by a list of index arrays
        """

        for i in range(len(self.__assoc)):
            self.__assoc[i] = self.__assoc[i][sorting[i]]

    def pick_first(self, sorting=None):
        """
        Returns the very first target associated with each fiber.
        """

        if sorting is None:
            sorting = len(self.__assoc) * [np.s_[:]]

        ix = np.array([ a[s][0] if a.size > 0 else -1 for a, s in zip(self.__assoc, sorting)], dtype=int)
        mask = (ix >= 0)
        return ix, mask

    def counts(self):
        """
        Returns the number of targets for each fiber
        """

        counts = np.array([ a.size for a in self.__assoc], dtype=int)
        return counts

    def histogram(self):
        """
        Calculates the histogram of target counts per fiber
        """
        hist = np.bincount(self.counts())
        return hist

    def plot_histogram(self, ax):        
        hist = self.histogram()
        ax.step(np.arange(hist.size) + 0.5, hist, where='mid')
        ax.set_xlabel("Number of stars per fiber")
        ax.set_ylabel("Number of fibers")
