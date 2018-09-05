import shutil
import os.path as path
from collections import namedtuple
from functools import wraps

from lsst.meas.algorithms import (LoadIndexedReferenceObjectsConfig, ScienceSourceSelectorConfig,
                                  ReferenceSourceSelectorConfig, DatasetConfig, IndexerRegistry,
                                  sourceSelectorRegistry)
from lsst.afw.table import unpackMatches
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase

__all__ = ["resolve_load_path", "generate_shard_paths", "get_shards", "toExposureF",
           "get_master_schema", "get_simple_catalog", "get_source_catalog",
           "persist_catalog", "get_exp_metadata", "get_exposure_id"]

################################################################################
############################  Path resolution  #################################
################################################################################
# Paths are a particular problem with this code that in the LSST Stack is handled 
# by the Butler. Without the Butler we are forced to manually create paths to 
# data we need.
def resolve_path_type(path):
    # looking at [:2] should cut out 'n' or : from s3n: or s3:
    # so don't use s3n because of replace will fail, generality?
    tmp = path.split("/")
    isS3 = tmp[0][:2] == "s3"
    bucket = tmp[2] if isS3 else None
    
    PathType = namedtuple("PathType", "isS3, bucket")
    return PathType(isS3, bucket)
            

def resolve_load_path(path):
    """Resolves whether provided path is an URI to an S3 bucket. If it is
    replaces the URI designation 's3://' with '~/s3-drive' where the bucket should
    be mounted before-hand. 

    Can't figure out how to read from bucket when it's not a Spark native type.
    """
    solvedPath = resolve_path_type(path)
    if solvedPath.isS3:
        path = path.replace("s3://"+solvedPath.bucket, "~/s3-drive")
    return path
#resolve_save_path = resolve_load_path


def generate_shard_paths(shardIds, refCatConf, refCatLoc):
    """Creates absolute paths to refcat shards from shard IDs. 

    Parameters
    ----------
    shardIds  : a list of shard IDs (integers)
    refCatConf : reference catalog configuration
    refCatLoc  : an absolute path to reference catalog (str)

    Returns
    ----------
    shard_locs : a list of paths to shards 
    """
    shard_locs = []
    ref_cat_path = resolve_load_path(refCatLoc)
    for shard_id in shardIds:
        rel_path  = path.join(refCatConf.ref_dataset_name, str(shard_id)+".fits")
        
        shard_loc = path.join(ref_cat_path, rel_path)
        shard_locs.append(shard_loc)

    return shard_locs


def generate_master_schema_path(refCatConf, refCatLoc):
    """Generate path to reference catalog master schema.

    Parameters
    ----------
    refCatConf : reference catalog configuration
    refCatLoc  : top level location where of the reference catalog.
    """
    rel_path = path.join(refCatConf.ref_dataset_name, "master_schema.fits")
    ref_cat_path = resolve_load_path(refCatLoc)
    master_schema_path = path.join(ref_cat_path, rel_path)
    return master_schema_path




################################################################################
##################  Reading files as LSST Stack objects  #######################
################################################################################
# Same goes when we convert things to objects the LSST Stack can understand. 
# Because that mapping was provided by the butler, which we don't have anymore, 
# we are forced to define the mappings ourselves. 

def toExposureF(exposure_loc):
    """Returns an exposure as a ExposureF object.
    
    Parameters
    ----------
    exposure_locs  : an absolute path to exposure
    """
    exposure_loc = resolve_load_path(exposure_loc)
    return afwImage.ExposureF(exposure_loc)


def get_shards(shard_locs): 
    """Creates a list of SourceCatalog objects from a list of paths to 
    reference catalog shards.

    Parameters
    ----------
    shard_locs  : a list of absolute paths to shards
    """
    shard_cats = []
    for shard_loc in shard_locs:
        shard_loc = resolve_load_path(shard_loc)
        cat = afwTable.SourceCatalog.readFits(shard_loc)
        shard_cats.append(cat)
    return shard_cats


def get_master_schema(master_schema_loc):
    """Read the reference catalog master schema as a Source Catalog.

    Parameters
    ----------
    master_schema_loc : absolute path to the master schema of reference catalog.
    """
    master_schema_loc = resolve_load_path(master_schema_loc)
    return afwTable.SourceCatalog.readFits(master_schema_loc)


def get_simple_catalog(cat_loc):
    """Read the fits as a SimpleCatalog.

    Parameters
    ----------
    cat_loc : absolute path the fits catalog
    """
    cat_loc = resolve_load_path(cat_loc)
    return afwTable.SimpleCatalog.readFits(cat_loc)


def get_source_catalog(cat_loc):
    """Read the fits as a SourceCatalog.

    Parameters
    ----------
    cat_loc : absolute path the fits catalog
    """
    cat_loc = resolve_load_path(cat_loc)
    return afwTable.SourceCatalog.readFits(cat_loc)



################################################################################
##############################  Persisting data.  ##############################
################################################################################
# Since we are not using the Butler any more there are no policies and naming conventions
# dictating where things should be kept and how exactly. For this reason we need to 
# wrap that functionality in somthing more Spark compliant.

def persist_catalog(location="", frmt=""):
    def persist(catWithId):
        cat, uniqueId = catWithId
        
        filename = frmt.format(uniqueId)
        saveloc = path.join(location, filename)
        saveloc = resolve_load_path(saveloc)
        
        tmp = path.join(path.expanduser("~"), filename)
        # cfitsio flips out if it tries to write to bucket
        cat.writeFits(tmp)
        # os.rename flips out when it tries to move between different file systems
        # shutil will work for python 3 but won't for older versions, shutil copy and 
        # os remove are the only general enough solution?
        shutil.move(tmp, path.expanduser(saveloc))
        return saveloc
    return persist




################################################################################
###############################  Metadata, Misc.  ##############################
################################################################################
# Since the butler is an essential part of the code for RefMatchTask we can't just use
# the original Task class which leaves a lot of methods (that read metadata or set-up 
# schemas, tables etc.) orphaned. 

def get_exp_metadata(exposure):
    """Return exposure's bounding box, wcs and, if existant, filter 
    and calibration data.    
    """
    expInfo = exposure.getInfo()    

    filterName = expInfo.getFilter().getName() or None
    if filterName == "_unknown_":
        filterName = None
        
    return pipeBase.Struct(
        bbox  = exposure.getBBox(),
        wcs   = expInfo.getWcs(),
        # calib objects not pickle-able, we don't seem to need it here but I wonder...
        #calib = expInfo.getCalib() if expInfo.hasCalib() else None,
        filterName = filterName,
    )


def get_exposure_id(exposure):
    """Returns exposure id.
    """
    expinf = exposure.getInfo()
    visinf = expinf.getVisitInfo()
    expid  = visinf.getExposureId()
    return expid


def add_flux_aliases(schema, refObjConf):
    """Add aliases for camera filter fluxes to the schema
    
    If self.config.defaultFilter then adds these aliases:
    camFlux:      <defaultFilter>_flux
    camFluxSigma: <defaultFilter>_fluxSigma, if the latter exists
    
    For each camFilter: refFilter in self.config.filterMap adds these aliases:
    <camFilter>_camFlux:      <refFilter>_flux
    <camFilter>_camFluxSigma: <refFilter>_fluxSigma, if the latter exists
    
    Parameters
    ----------
    schema     : schema to which flux aliases will be appended
    refObjConf : a LoadReferenceObject configuration that defines the pairing
                 ("translation") between reference catalog and science catalog
                 schemas.
    """
    aliasMap = schema.getAliasMap()
    
    def add_filter_alias(filterName, refFilterName):
        """Add aliases for a single filter
        
        Parameters
        ----------
        filterName    :  camera filter name, or "". Defauts to  <filterName>_camFlux 
                         or camFlux if filterName is None
        refFilterName :  reference filter name; <refFilterName>_flux must exist
        """
        camFluxName = filterName + "_camFlux" if filterName is not None else "camFlux"
        refFluxName = refFilterName + "_flux"
        if refFluxName not in schema:
            raise RuntimeError("Unknown reference filter %s" % (refFluxName,))
        aliasMap.set(camFluxName, refFluxName)
        refFluxErrName = refFluxName + "Sigma"
        if refFluxErrName in schema:
            camFluxErrName = camFluxName + "Sigma"
            aliasMap.set(camFluxErrName, refFluxErrName)

    if refObjConf.defaultFilter:
        add_filter_alias(None, refObjConf.defaultFilter)

    for filterName, refFilterName in refObjConf.filterMap.items():
        add_filter_alias(filterName, refFilterName)


def get_ref_flux_field(schema, filterName=None):
    """Get name of flux field in schema.

    if <filterName> is specified:
        return <filterName>_camFlux if present
        else return <filterName>_flux if present (camera filter name matches reference filter name)
        else throw RuntimeError
    else:
        return camFlux, if present,
        else throw RuntimeError

    Parameters
    ----------
    schema     : reference catalog schema
    filterName : name of camera filter
    """
    if not isinstance(schema, afwTable.Schema):
        raise RuntimeError("schema=%s is not a schema" % (schema,))
    if filterName:
        fluxFieldList = [filterName + "_camFlux", filterName + "_flux"]
    else:
        fluxFieldList = ["camFlux"]
    for fluxField in fluxFieldList:
        if fluxField in schema:
            return fluxField

    raise RuntimeError("Could not find flux field(s) %s" % (", ".join(fluxFieldList)))


def add_centroid_columns(refCat):
    """Adds columns required for matching to the reference catalog and returns
    the expanded catalog. 

    Ads columns:
        - centroid_x
        - centroid_y
        - hasCentroid

    Parameters
    ----------
    refCat : a reference catalog to which we want to append columns
    """
    # add and initialize centroid and hasCentroid fields (these are added
    # after loading to avoid wasting space in the saved catalogs)
    # the new fields are automatically initialized to (nan, nan) and False
    mapper = afwTable.SchemaMapper(refCat.schema, True)
    mapper.addMinimalSchema(refCat.schema, True)
    mapper.editOutputSchema().addField("centroid_x", type=float)
    mapper.editOutputSchema().addField("centroid_y", type=float)
    mapper.editOutputSchema().addField("hasCentroid", type="Flag")

    expandedCat = afwTable.SimpleCatalog(mapper.getOutputSchema())
    expandedCat.extend(refCat, mapper=mapper)
    
    return expandedCat

