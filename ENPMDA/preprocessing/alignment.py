"""\
=========
Alignment
=========
The `Alignment` module contains a template for running alignment
on a `Universe`.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align
from loguru import logger

class AlignmentBase:
    """
    Base class for alignment of a `Universe`. Run the `process_universe` method
    to align the `Universe`.

    Properties
    ----------
    universe : MDAnalysis.Universe
        The Universe to align.
    output_prefix : str
        The prefix for the output file.
    """
    def __init__(self):
        self._universe = None
        self._output_prefix = None

    @property
    def universe(self):
        return self._universe

    @universe.setter
    def universe(self, value):
        if not isinstance(value, mda.Universe):
            raise TypeError("The universe must be an MDAnalysis Universe.")
        self._universe = value
    
    @property
    def output_prefix(self):
        return self._output_prefix

    @output_prefix.setter
    def output_prefix(self, value):
        if not isinstance(value, str):
            raise TypeError("The output prefix must be a string.")
        self._output_prefix = value

    def process_universe(self):
        raise NotImplementedError("This method must be implemented in a subclass.")
    
    # create a copy method with the same class instance but empty universe and output_prefix
    def copy(self):
        """
        Return a shallow copy of this instance with universe and output_prefix reset.
        Works even if subclass has a different __init__.
        """
        # Create a new empty instance of the same class without calling __init__
        new_obj = self.__class__.__new__(self.__class__)

        # Copy all attributes shallowly
        new_obj.__dict__ = self.__dict__.copy()

        # Reset universe and output_prefix
        new_obj._universe = None
        new_obj._output_prefix = None

        return new_obj


class CalphaAlignment(AlignmentBase):
    def process_universe(self):
        aligner = align.AlignTraj(
                    mobile=self.universe,
                    reference=self.universe,
                    select="protein and name CA",
                    filename=self.output_prefix + ".xtc")
        aligner.run()
        logger.info("C-alpha alignment complete.")