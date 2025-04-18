-- USE HLSP_PS1_PSC

DECLARE @ra AS float = 227.27229069;
DECLARE @dec AS float = 67.23449905;
DECLARE @rad AS float = 120; 

SELECT 
    p.objID as obj_id,
    p.raMean AS RA, p.decMean AS Dec,
    gPSFMag, gPSFMagErr, rPSFMag, rPSFMagErr, iPSFMag, iPSFMagErr, zPSFMag, zPSFMagErr, yPSFMag, yPSFMagErr
INTO Mydb.PS1_GAIA_point_source_Ursaminor_1
FROM pointsource_magnitudes_view AS p
INNER JOIN fGetNearbyObjEq(@ra, @dec, @rad) AS r ON p.objid=r.objid
WHERE p.ps_score > 0.9 AND primaryDetection = 1
    AND gPSFMag != -999 AND rPSFMag != -999 AND zPSFMag != -999
       
CREATE CLUSTERED INDEX IX_PS1_GAIA_point_source_Ursaminor_1 ON PS1_GAIA_point_source_Ursaminor_1 (obj_id)