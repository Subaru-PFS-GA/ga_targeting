from .ursaminor import UrsaMinor as __UrsaMinor
from .sculptor import Sculptor as __Sculptor
from .draco import Draco as __Draco
from .fornax import Fornax as __Fornax
from .bootes import Bootes as __Bootes
from .ngc6822 import NGC6822 as __NGC6822

# Singleton instances

UrsaMinor = __UrsaMinor()
Sculptor = __Sculptor()
Draco = __Draco()
Fornax = __Fornax()
Bootes = __Bootes()
NGC6822 = __NGC6822()

GALAXIES = {
    'umi': UrsaMinor,
    'scl': Sculptor,
    'dra': Draco,
    'for': Fornax,
    'booi': Bootes,
    'ngc6822': NGC6822,
}
