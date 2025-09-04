# Prefixes used to prefix the target index used internally to identify
# science targets. These are chosen to be large enough to not overlap with
# the target indices. Note that these are not targetids or objid used in
# source databases, those number are always left intact.

ID_PREFIX_BOOTES    = 0x0100000000
ID_PREFIX_DRACO     = 0x0200000000
ID_PREFIX_FORNAX    = 0x0300000000
ID_PREFIX_NGC6822   = 0x0400000000
ID_PREFIX_SCULPTOR  = 0x0500000000
ID_PREFIX_URSAMINOR = 0x0600000000

# Outer disk fields
ID_PREFIX_OD_L90_B28   = 0x1000000000
ID_PREFIX_OD_L90_B29   = 0x2000000000
ID_PREFIX_OD_L90_BM28  = 0x3000000000
ID_PREFIX_OD_L90_BM29  = 0x4000000000
ID_PREFIX_OD_L90_B25   = 0x5000000000
ID_PREFIX_OD_L90_B27   = 0x6000000000
ID_PREFIX_OD_L90_B16   = 0x7000000000

# Cross-calibration fields
ID_PREFIX_CC_RA288_DECM11 = 0x10000000000
ID_PREFIX_CC_RA288_DECM17 = 0x20000000000
ID_PREFIX_CC_RA288_DECM22 = 0x30000000000
ID_PREFIX_CC_RA336_DECM12 = 0x40000000000


# M31 fields, sectors will have four bits following
# the 1 set
ID_PREFIX_M31       = 0x100000000000000
ID_PREFIX_M31_E0    = ID_PREFIX_M31 | (0x0001 << (32))
ID_PREFIX_M31_W0    = ID_PREFIX_M31 | (0x0002 << (32))
ID_PREFIX_M31_GSS0  = ID_PREFIX_M31 | (0x0003 << (32))
ID_PREFIX_M31_NWS0  = ID_PREFIX_M31 | (0x0004 << (32))