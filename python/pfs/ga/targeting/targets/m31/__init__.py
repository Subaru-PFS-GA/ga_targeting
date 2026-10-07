from .m31 import M31
from .m33 import M33

M31_SECTORS = {
    'all': M31(),
    'm33': M33(),
}

for sector, field, _, _, _, _, _ in M31.FIELDS:
    if sector not in M31_SECTORS:
        M31_SECTORS[sector] = M31(sector=sector)
    if field not in M31_SECTORS:
        M31_SECTORS[field] = M31(sector=sector, field=field)
        
for sector, field, _, _, _, _, _ in M33.FIELDS:
    if sector not in M31_SECTORS:   
        M31_SECTORS[sector] = M33(sector=sector)
    if field not in M31_SECTORS:
        M31_SECTORS[field] = M33(sector=sector, field=field)