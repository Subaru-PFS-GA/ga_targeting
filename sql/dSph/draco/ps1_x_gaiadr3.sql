-- USE GAIA_DR3

SELECT
    m.obj_id, g.source_id, ref_epoch,
    g.ra, g.dec,
    g.parallax, g.parallax_error,
    g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
    g.radial_velocity, g.radial_velocity_error,
    n.angular_distance,
    gPSFMag, gPSFMagErr, rPSFMag, rPSFMagErr, iPSFMag, iPSFMagErr, zPSFMag, zPSFMagErr, yPSFMag, yPSFMagErr,
    phot_g_mean_flux, phot_g_mean_flux_error,
    phot_bp_mean_flux, phot_bp_mean_flux_error,
    phot_rp_mean_flux, phot_rp_mean_flux_error,
    number_of_neighbours, number_of_mates
INTO PS1_GAIA_point_source_Draco_2
FROM MYDB.PS1_GAIA_point_source_Draco_1 m 
INNER JOIN panstarrs1_best_neighbour n
    ON m.obj_id = n.clean_panstarrs1_oid
INNER JOIN gaia_source g
    ON g.source_id = n.source_id