from .ursaminor import UrsaMinor as __UrsaMinor
from .sculptor import Sculptor as __Sculptor
from .draco import Draco as __Draco
from .fornax import Fornax as __Fornax
from .bootes import Bootes as __Bootes
from .ngc6822 import NGC6822 as __NGC6822
from .sextans import Sextans as __Sextans

# Singleton instances

UrsaMinor = __UrsaMinor()
Sculptor = __Sculptor()
Draco = __Draco()
Fornax = __Fornax()
Bootes = __Bootes()
NGC6822 = __NGC6822()
Sextans = __Sextans()

GALAXIES = {
    'umi': UrsaMinor,
    'ursaminor': UrsaMinor,
    'scl': Sculptor,
    'sculptor': Sculptor,
    'dra': Draco,
    'draco': Draco,
    'for': Fornax,
    'fornax': Fornax,
    'booi': Bootes,
    'bootes': Bootes,
    'bootesi': Bootes,
    'ngc6822': NGC6822,
    'sex': Sextans,
    'sextans': Sextans,
}
