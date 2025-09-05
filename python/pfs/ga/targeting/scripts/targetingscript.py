import os
from collections.abc import Iterable
import pandas as pd
import numpy as np
import json

from pfs.datamodel import TargetType

from pfs.ga.common.io import ObservationSerializer
from pfs.ga.common.photometry import PhotometryEncoder, PhotometryDecoder

from ..config.netflow import NetflowConfig
from ..instrument import SubaruPFI
from ..targets.dsph import GALAXIES as DSPH_FIELDS
from ..targets.m31 import M31_SECTORS
from . import Script

from ..setup_logger import logger

class TargetingScript(Script):
    """
    Base class that implements common function for targeting scripts
    """

    def __init__(self):
        super().__init__()

        self._indir = None
        self._outdir = None
        self._field = None
        self._config = None

    def __get_config(self):
        return self._config
    
    config = property(__get_config)

    def _add_args(self):
        super()._add_args()

    def _create_field_from_args(self, args):
        if self.is_arg('dsph', args):
            self._field = DSPH_FIELDS[self.get_arg('dsph', args)]

        if self.is_arg('m31', args):
            self._field = M31_SECTORS[self.get_arg('m31', args)]

    def _create_config_from_field(self):
        # If a field is specified, load its default configuration      
        if self._field is not None:
            self._config = self._field.get_netflow_config()
        else:
            self._config = NetflowConfig.default()

    def _load_config_files(self, args):
        # Load the configuration template files and merge with the default config
        config_files = self.get_arg('config', args)
        self._config.load(config_files, ignore_collisions=True)
        logger.info(f'Loaded {len(config_files)} config files from {config_files}.')
        logger.info(f'Found {len(self._config.targets)} target lists in the configuration.')

    def _get_output_config_path(self):
        command = self.get_command_name()
        path = os.path.join(self._outdir, f'{command}_{self.timestamp}.config')
        return path

    def _dump_settings(self):
        super()._dump_settings()

        # Save the active configuration to the output directory
        path = self._get_output_config_path()
        self._config.save(path)

        logger.debug(f'Configuration saved to `{os.path.abspath(path)}`.')

    def _create_instrument(self):
        return SubaruPFI(instrument_options=self._config.instrument_options)
    
    def _generate_pointings(self, stage=None):
        # The list of pointings can be defined in the config file, which then
        # override the pointings defined for the field in the library source code.

        if stage is not None:
            if not isinstance(stage, Iterable):
                stage = [ stage ]
            stage = set(stage)

        # Get the pointing for the config file or from the field, if defined
        if self._config.pointings is not None:
            pp = [ p.get_pointing() for p in self._config.pointings ]
        elif self._field is not None:
            # Load from the class
            pp = self._field.get_pointings(SubaruPFI)
        else:
            raise NotImplementedError()
        
        # The number of visits can be overridden from the command-line argument,
        # see __init_from_args()
        nvisits = self._config.field.nvisits
        nrepeats = self._config.field.nrepeats
        exp_time = self._config.field.exp_time
        obs_time = self._config.field.obs_time
        
        # Create new pointing objects with obs_time, exp_time etc.
        pointings = []
        for p in pp:   
            if stage is None or p.stage in stage:
                if obs_time is not None:
                    p.obs_time = obs_time
                if exp_time is not None:
                    p.exp_time = exp_time
                if nvisits is not None:
                    p.nvisits = nvisits
                if nrepeats is not None:
                    p.nrepeats = nrepeats
                pointings.append(p)

        if len(pointings) == 0:
            raise ValueError(f'No pointings found, cannot continue.')
        else:
            logger.info(f'Found {len(pointings)} pointings for observation stage {stage}. '
                        f'Number of visits is {nvisits} with exposure time {exp_time} s. '
                        f'Every visit will be repeated {nrepeats} times. '
                        f'Observation time is {obs_time} UTC.')
            
        return pointings

    @staticmethod
    def load_target_list(key, target_list_config, fn=None):
        """
        Loads an input target list and performs column mapping as well as the discovery of the
        available photometric systems.
        """
        
        if target_list_config.reader is not None:
            # TODO: implement custom reader function
            raise NotImplementedError()
        
        fn = fn if fn is not None else os.path.expandvars(target_list_config.path)
        
        logger.info(f'Reading target list from `{fn}`.')

        reader = ObservationSerializer(
            catalog_name=key,
            columns = target_list_config.columns,
            column_map = target_list_config.column_map,
            data_types = target_list_config.data_types,
            value_map = target_list_config.value_map,
            index = target_list_config.index,
            filters = target_list_config.photometry.filters if target_list_config.photometry is not None else None,
            bands = target_list_config.photometry.bands if target_list_config.photometry is not None else None,
            limits = target_list_config.photometry.limits if target_list_config.photometry is not None else None,
            mask = target_list_config.mask,
            kwargs = target_list_config.reader_args,
        )
        catalog = reader.read(fn)
        catalog.name = key
        
        logger.info(f'Read {catalog.shape[0]} targets from target list.')

        return catalog


    def _get_id_prefix(self):
        if self._config.field.id_prefix is not None:
            id_prefix = self._config.field.id_prefix
        elif self._field is not None and self._field.id_prefix is not None:
            id_prefix = self._field.id_prefix
        else:
            id_prefix = 0

        return id_prefix

    def _get_preprocessed_target_list_path(self, key, dir):
        return os.path.join(dir, f'{self._config.field.key}_targets_{key}.feather')

    def _save_preprocessed_target_lists(self, target_lists, dir):
        for k, target_list_config in self._config.targets.items():
            path = self._get_preprocessed_target_list_path(k, dir)
            self.__save_preprocessed_target_list(k, target_lists[k], path)

    def __save_preprocessed_target_list(self, key, target_list, path):
        logger.info(f'Saving preprocessed target list `{key}` to `{path}`.')
        s = ObservationSerializer()
        s.write(target_list, path)

        fn, ext = os.path.splitext(path)
        path = f'{fn}.json'
        logger.info(f'Saving photometry info for target list `{key}` to `{path}`.')
        with open(path, 'w') as f:
            json.dump(target_list.photometry, f,
                      indent=4, cls=PhotometryEncoder)

    def _load_preprocessed_target_list(self, key, path):
        """
        Load a preprocessed target list. This target list already has column names
        mapped to internal column names and the flux columns unwrapped so this function
        is a simplified version of `load_target_list`.
        """

        logger.info(f'Loading preprocessed target list `{key}` from `{path}`.')
        s = ObservationSerializer()
        target_list = s.read(path)
        target_list.name = key

        fn, ext = os.path.splitext(path)
        path = f'{fn}.json'
        logger.info(f'Loading photometry info for target list `{key}` from `{path}`.')
        with open(path, 'r') as f:
            phot = json.load(f, cls=PhotometryDecoder)
            target_list._set_photometry(phot)

        return target_list
    
    def _get_design_list_path(self):
        raise NotImplementedError()

    def _save_design_list(self, netflow, designs):
        df = pd.DataFrame({
            'stage': [ v.pointing.stage for v in netflow.visits ],
            'pointing_idx': [ v.pointing_idx for v in netflow.visits ],
            'visit_idx': [ v.visit_idx for v in netflow.visits ],
            'priority': [ v.pointing.priority for v in netflow.visits ],
            'nrepeats': [ v.pointing.nrepeats for v in netflow.visits ],
            'ra': [ v.pointing.ra for v in netflow.visits ],
            'dec': [ v.pointing.dec for v in netflow.visits ],
            'posang': [ v.pointing.posang for v in netflow.visits ],
            'obs_time': [ v.pointing.obs_time.iso for v in netflow.visits ],
            'exp_time': [ v.pointing.exp_time.value for v in netflow.visits ],
            'pfsDesignId': [ d.pfsDesignId for d in designs]
        })

        fn = self._get_design_list_path()
        df.to_feather(fn)

        logger.info(f'Saved list of designs to `{fn}`.')

    def _get_fiber_assignments_all_path(self, dir):
        return os.path.join(dir, f'{self._config.field.key}_assignments_all.feather')

    def _save_fiber_assignments_all(self, fiber_assignments_all):
        # This dataframe contains columns with dtype `object`` that contain the list of
        # magnitudes, fluxes, filter names, etc.
        
        path = self._get_fiber_assignments_all_path(self._outdir)
        logger.info(f'Saving assignments joined with the input catalogs to `{path}`.')
        fiber_assignments_all.to_feather(path)

    def _load_fiber_assignments_all(self, dir):
        path = self._get_fiber_assignments_all_path(dir)
        logger.info(f'Loading assignments joined with the input catalogs from `{path}`.')
        return pd.read_feather(path)

    def _create_pfsDesign_visit(self, visit, assignments,
                               design_name='',
                               arms='bmn'):
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

        # Use targetid for sky and cal targets and use target_idx for sci targets
        target_type = np.array(fiber_assignments['target_type'].fillna(-1).astype(np.int32))
        sci_mask = target_type == TargetType.SCIENCE
        target_idx = np.array(fiber_assignments['__target_idx'].fillna(-1).astype(np.int64))
        obj_id = np.array(fiber_assignments['targetid'].fillna(-1).astype(np.int64))
        
        id_prefix = self._get_id_prefix()
        obj_id[sci_mask] = (id_prefix | target_idx[sci_mask])

        kwargs = dict(
            designName = design_name,
            variant = 0,          # not used by makePfsDesign
            designId0 = 0,        # not used by makePfsDesign

            raBoresight = visit.pointing.ra,
            decBoresight = visit.pointing.dec,
            posAng = visit.pointing.posang,
            arms = ''.join(arms),
            fiberId = np.array(fiber_assignments.index.astype(np.int32)),
            tract = np.array(fiber_assignments['tract'].fillna(-1).astype(np.int32)),
            patch = np.array(fiber_assignments['patch'].fillna('0,0').astype(str)),
            ra = np.array(fiber_assignments['RA'].astype(np.float64)),
            dec = np.array(fiber_assignments['Dec'].astype(np.float64)),
            catId = np.array(fiber_assignments['catid'].fillna(-1).astype(np.int32)),
            objId = obj_id,
            targetType = target_type,
            fiberStatus = np.array(fiber_assignments['fiber_status'].fillna(-1).astype(np.int32)),
            epoch = np.array(fiber_assignments['epoch'].astype(str)),
            pmRa = np.array(fiber_assignments['pmra'].astype(np.float64)),
            pmDec = np.array(fiber_assignments['pmdec'].astype(np.float64)),
            parallax = np.array(fiber_assignments['parallax'].astype(np.float64)),
            proposalId = np.array(fiber_assignments['proposalid'].astype(str)),
            
            obCode = np.array(fiber_assignments['obcode']),    # TODO: this should be unique for each design

            pfiNominal = np.stack([
                 fiber_assignments['fp_x'].astype(float),
                 fiber_assignments['fp_y'].astype(float)], axis=-1),

            guideStars = None,
        )

        # Substitute na values with empty lists, it's a bit hacky because pandas
        # does not support passing lists to fillna
        for c in ['filter', 'fiber_flux', 'fiber_flux_err', 'psf_flux', 'psf_flux_err', 'total_flux', 'total_flux_err']:
            isna = fiber_assignments[c].isna()
            fiber_assignments.loc[fiber_assignments.index[isna], c] = fiber_assignments.loc[fiber_assignments.index[isna], c].apply(lambda x: [])

        # Add the fluxes
        # Convert the rows of column `filter` into a list
        kwargs['filterNames'] = list(fiber_assignments['filter'])
        for prefix in ['fiber', 'psf', 'total']:
                kwargs[f'{prefix}Flux'] = list(fiber_assignments[f'{prefix}_flux'])
                kwargs[f'{prefix}FluxErr'] = list(fiber_assignments[f'{prefix}_flux_err'])

        # Calculate the design ID hash from the fibers and coordinates
        kwargs['pfsDesignId'] = calculate_pfsDesignId(kwargs['fiberId'], kwargs['ra'], kwargs['dec'])

        pfsDesign = PfsDesign(**kwargs)

        return pfsDesign