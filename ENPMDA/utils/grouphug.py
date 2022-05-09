# From https://groups.google.com/g/mdnalysis-discussion/c/umDpvbCmQiE/m/FKtNClazAwAJ
# Author: Richard Gowers

from MDAnalysis.transformations.base import TransformationBase
import numpy as np


class GroupHug(TransformationBase):
    def __init__(self, center, *others):
        super().__init__(max_threads=1,
                         parallelizable=True)
        self.c = center
        self.o = others

    @staticmethod
    def calc_restoring_vec(ag1, ag2):
        box = ag1.dimensions[:3]
        dist = ag1.center_of_mass() - ag2.center_of_mass()

        return box * np.rint(dist / box)

    def _transform(self, ts):
        # loop over other atomgroups and shunt them into nearest image to
        # center
        for i in self.o:
            rvec = self.calc_restoring_vec(self.c, i)
            i.translate(+rvec)
        return ts
