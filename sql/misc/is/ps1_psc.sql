-- USE HLSP_PS1_PSC

DECLARE @ra AS float = 150.05;
DECLARE @dec AS float = 2.2;
DECLARE @rad AS float = 100; 

SELECT 
    p.objID as obj_id,
    p.raMean AS RA, p.decMean AS Dec,
    gPSFMag, gPSFMagErr, rPSFMag, rPSFMagErr, iPSFMag, iPSFMagErr, zPSFMag, zPSFMagErr, yPSFMag, yPSFMagErr
INTO Mydb.PS1_GAIA_point_source_ra150dec2
FROM pointsource_magnitudes_view AS p
INNER JOIN fGetNearbyObjEq(@ra, @dec, @rad) AS r ON p.objid=r.objid
WHERE p.ps_score > 0.9 AND primaryDetection = 1
    AND gPSFMag != -999 AND rPSFMag != -999 AND zPSFMag != -999
       