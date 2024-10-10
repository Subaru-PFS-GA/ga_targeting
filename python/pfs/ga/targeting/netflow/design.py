import numpy as np

from collections.abc import Iterable

class Design():
    """
    Utility class to create PfsDesign object from a list of targets.
    """

    def join_catalogs(assignments, catalogs):
        catalogs = [catalogs] if not isinstance(catalogs, Iterable) else catalogs

        for catalog in catalogs:
            assignments = assignments.merge(catalog, on='targetid', how='left')

    def get_pfsDesign_visit(visit, assignments):
        """
        Generate a PfsDesign object for a given visit.
        """

        # TODO: add proposal_id, obCode postfix or prefix, designName, variant, designId0
        #       what about the various fluxes? we have PSF flux only
        #       tract, patch from coordinates
        #       what about guide stars?

        from pfs.datamodel import PfsDesign
        from pfs.datamodel.utils import calculate_pfsDesignId
        
        # Filter down assignment list to the current visit
        mask = (assignments['visit_idx'] == visit.visit_idx) & \
               (assignments['pointing_idx'] == visit.pointing_idx)
        
        fiber_assignments = assignments[mask].set_index(['fiberid'])
        fiber_assignments.sort_index(inplace=True)

        nfibers = len(fiber_assignments)
                
        kwargs = dict(
            designName = 'ga_{galaxy.ID}',
            variant = 0,
            designId0 = 0,

            raBoresight = visit.pointing.ra,
            decBoresight = visit.pointing.dec,
            posAng = visit.pointing.posang,
            arms = 'bmn',
            fiberId = np.array(fiber_assignments.index, dtype=np.int32),
            tract = np.array(fiber_assignments['tract'], dtype=np.int32),
            patch = np.array(fiber_assignments['patch']),
            ra = np.array(fiber_assignments['RA'], dtype=np.float64),
            dec = np.array(fiber_assignments['Dec'], dtype=np.float64),
            catId = np.array(fiber_assignments['catid'], dtype=np.int32),
            objId = np.array(fiber_assignments['targetid'], dtype=np.int64),
            targetType = np.array(fiber_assignments['target_type'], dtype=np.int64),
            fiberStatus = np.array(fiber_assignments['fiber_status'], dtype=np.int64),
            epoch = np.array(fiber_assignments['epoch']),
            pmRa = np.array(fiber_assignments['pmra'], dtype=np.float64),
            pmDec = np.array(fiber_assignments['pmdec'], dtype=np.float64),
            parallax = np.array(fiber_assignments['parallax'], dtype=np.float64),
            proposalId = np.array(fiber_assignments['proposalid']),
            
            obCode = np.array(fiber_assignments['obcode']),    # TODO: this should be unique for each design

            pfiNominal = np.stack([ fiber_assignments['fp_x'],  fiber_assignments['fp_y']], axis=-1).astype(float),

            guideStars = None,
        )

        # Add the fluxes
        # Convert the rows of column `filter` into a list
        kwargs['filterNames'] = list(fiber_assignments['filter'])
        for prefix in ['fiber', 'psf', 'total']:
                kwargs[f'{prefix}Flux'] = list(fiber_assignments[f'{prefix}_flux'])
                kwargs[f'{prefix}FluxErr'] = list(fiber_assignments[f'{prefix}_flux_err'])

        # Calculate the design ID hash from the fibers and coordinates
        kwargs['pfsDesignId'] = calculate_pfsDesignId(kwargs['fiberId'], kwargs['ra'], kwargs['dec'])

        return PfsDesign(**kwargs)
    
    def get_pfsDesign_all(self, filters, assignments=None):
        """
        Generate a list of PfsDesign objects for all visits.
        """

        designs = []
        for visit in self.__visits:
            designs.append(self.get_pfsDesign_visit(visit, filters, assignments=assignments))

        return designs