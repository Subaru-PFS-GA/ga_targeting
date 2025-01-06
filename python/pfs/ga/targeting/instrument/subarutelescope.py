from pfs.utils.coordinates import DistortionCoefficients as DCoeff

from .instrument import Instrument

class SubaruTelescope(Instrument):
    """
    Implements functions specific to the Subaru Telescope. Most of these
    are just wrappers around the PFS libraries.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
