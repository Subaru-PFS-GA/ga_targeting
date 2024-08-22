from .ursaminor import UrsaMinor as __UrsaMinor
from .sculptor import Sculptor as __Sculptor

# Singleton instances

UrsaMinor = __UrsaMinor()
Sculptor = __Sculptor()

GALAXIES = {
    'umi': UrsaMinor,
    'scl': Sculptor
}
