class CobraAngleFlags:
    SOLUTION_OK           = 0x0001  # 1 if the solution is valid
    IN_OVERLAPPING_REGION = 0x0002  # 1 if the position in overlapping region
    PHI_NEGATIVE          = 0x0004  # 1 if phi angle is negative(phi CCW limit < 0)
    PHI_BEYOND_PI         = 0x0008  # 1 if phi angle is beyond PI(phi CW limit > PI)
    TOO_CLOSE_TO_CENTER   = 0x0010  # 1 if the position is too close to the center
    TOO_FAR_FROM_CENTER   = 0x0020  # 1 if the position is too far from the center