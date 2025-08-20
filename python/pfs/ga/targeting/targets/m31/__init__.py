from .m31 import M31

M31_SECTORS = {
    'all': M31(),
}

for sector, field, _, _, _, _, _ in M31.FIELDS:
    if sector not in M31_SECTORS:
        M31_SECTORS[sector] = M31(sector=sector)
    if field not in M31_SECTORS:
        M31_SECTORS[field] = M31(sector=sector, field=field)