import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import IndexerRegistry, sourceSelectorRegistry
from lsst.meas.astrom import MatchPessimisticBTask
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

from . import utils

__all__ = ["create_reference_catalogs", "match_sources_with_metadata",
           "match_sources", "resolve_shard_ids"]

################################################################################
##############################  ID Resolvers  ##################################
################################################################################

def calculate_circle(bbox, wcs, pixelMargin):
    """Computes on-sky center and radius of search region

    Parameters
    ----------
    bbox :  bounding box for pixels (an lsst.geom.Box2I or Box2D)
    wcs  :  WCS (an lsst.afw.geom.SkyWcs)
    pixelMargin : padding in pixels by which the bbox will be expanded

    Returns
    ----------
    coord  : ICRS center of the search region (lsst.geom.SpherePoint)
    radius : the radius of the search region (lsst.geom.Angle)
    bbox   : the bounding box used to compute the circle (lsst.geom.Box2D)
    """
    bbox = afwGeom.Box2D(bbox) 
    bbox.grow(pixelMargin)
    coord  = wcs.pixelToSky(bbox.getCenter())
    radius = max(coord.separation(wcs.pixelToSky(pp)) for pp in bbox.getCorners())
    return pipeBase.Struct(coord=coord, radius=radius, bbox=bbox)


def resolve_circle2shard_ids(refCatConf, circle):
    """Resolves IDs of shards overlapping an on-sky circular region. 

    Parameters
    ----------
    refCatConf : configuration of the reference catalog 
                 (an lsst.meas.algorithms.DatasetConfig object)
    circle     : an lsst.pipe.base.Struct object containing ICRS center of
                 search region (an lsst.geom.SpherePoint) and radius
                 of the search region (an lsst.geom.Angle) 

    Returns
    ----------
    shard_ids     : a list of integer IDs of reference catalog shards 
                    overlapping the region
    boundary_mask : a boolean array indicating whether the shard touches the 
                    boundary (True) or is fully contained (False)
    """
    ref_dataset_name = refCatConf.ref_dataset_name
    indexer = IndexerRegistry[refCatConf.indexer.name](refCatConf.indexer.active)

    id_list, boundary_mask = indexer.get_pixel_ids(circle.coord, circle.radius)

    shard_ids = []
    for pixel_id in id_list:
        # how would you ask if datasetExists without butler?
        shard_ids.append(pixel_id)
        
    return shard_ids, boundary_mask


def resolve_bbox2shard_ids(refCatConf, bbox, wcs, pixelMargin=300):
    """Resolves IDs of shards overlapping an on-sky bounding box.

    Parameters
    ----------
    refCatConf : configuration of the reference catalog 
                 (an lsst.meas.algorithms.DatasetConfig object)
    bbox       : bounding box of the region of interest (lsst.geom.Box2D object)
    wcs        : WCS defining the bbox coordinate system (lsst.afw.geom.SkyWcs)
    pixelMargin: padding to add to 4 all edges of the bounding box (pixels)
                 default 300

    Returns
    ----------
    circle    : an lsst.pipe.base.Struct object containing ICRS center of
                search region (an lsst.geom.SpherePoint) and radius
                of the search region (an lsst.geom.Angle) 
    shard_ids : a list of integer IDs of reference catalog shards overlapping
                the region
    boundary_mask : a boolean array indicating whether the shard touches the 
                    boundary (True) or is fully contained (False)
    """
    circle = calculate_circle(bbox, wcs, pixelMargin)
    shard_ids, boundary_mask = resolve_circle2shard_ids(refCatConf, circle)
    return circle, shard_ids, boundary_mask


def resolve_shard_ids(refCatConf, exposure, **kwargs):
    """Resolves IDs of shards overlapping an exposure.

    Parameters
    ----------
    refCatConf : configuration of the reference catalog 
                 (an lsst.meas.algorithms.DatasetConfig object)
    exposure   : an exposure of the sky (lsst.afw.ExposureF)
    pixelMargin: padding to add to 4 all edges of the bounding box (pixels)
                 default 300
    Returns
    ----------
    circle    : an lsst.pipe.base.Struct object containing ICRS center of
                search region (an lsst.geom.SpherePoint) and radius
                of the search region (an lsst.geom.Angle) 
    shard_ids : a list of integer IDs of reference catalog shards overlapping
                the region
    boundary_mask : a boolean array indicating whether the shard touches the 
                    boundary (True) or is fully contained (False)
    """

    meta = utils.get_exp_metadata(exposure)
    return resolve_bbox2shard_ids(refCatConf, bbox=meta.bbox, 
                                  wcs=meta.wcs, **kwargs)






################################################################################
#######################  Reference Catalog Creation  ###########################
################################################################################

def trim2circle(catalogShard, ctrCoord, radius):
    """Trim a catalog to a circular aperture.
    
    Parameters
    ----------
    catalog_shard : SourceCatalog to be trimmed
    ctrCoord      : afw.Coord to compare each record to
    radius        : afwGeom.Angle indicating maximume separation
    """
    temp_cat = type(catalogShard)(catalogShard.schema)
    for record in catalogShard:
        if record.getCoord().separation(ctrCoord) < radius:
            temp_cat.append(record)
    return temp_cat


def trim2bbox(refCat, bbox, wcs):
    """Remove objects outside a given pixel-based bbox and set centroid and
    hasCentroid fields.
    
    Parameters
    ----------
    refCat : a catalog of objects (an lsst.afw.table.SimpleCatalog, or other
             table type that has fields "coord", "centroid" and "hasCentroid").
    bbox   : pixel region (an afwImage.Box2D)
    wcs    : WCS used to convert sky position to pixel position (an lsst.afw.math.WCS)
    """
    afwTable.updateRefCentroids(wcs, refCat)
    centroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
    retStarCat  = type(refCat)(refCat.table)
    for star in refCat:
        point = star.get(centroidKey)
        if bbox.contains(point):
            retStarCat.append(star)
    return retStarCat


def merge_and_mask2circle(circle, shards, boundaryMask, refCatConf,
                          refCatLoc, refObjConf):
    """ Reads the master schema of the reference catalog and fills it with
    sources that intersect the desired sky-circle for each shard in shards list. 
    Returned catalog is not neccessarily contiguous.

    Parameters
    ----------
    circle       : an lsst.pipe.base.Struct object containing ICRS center of
                   search region (an lsst.geom.SpherePoint) and radius
                   of the search region (an lsst.geom.Angle) 
    shards       : a list of shards (SourceCataog objects)
    boundaryMask : a boolean array indicating whether the shard touches the 
                    boundary (True) or is fully contained (False)
    refCatConf   : reference catalog configuration
    location     : top level location of the reference catalog
    """
    master_schema_loc = utils.generate_master_schema_path(refCatConf, refCatLoc)
    refCat = utils.get_master_schema(master_schema_loc)
    utils.add_flux_aliases(refCat.schema, refObjConf)

    # flux field identifies the filter of the exposure on which the trims were
    # made because matcher needs to know which columns to compare. Removing it 
    # means it needs to be specified manually during the matching.
    # The "problem" being that exposures come with the information in which filter 
    # they were made and remove any ambiguity when matching and therefore should
    # not be the users responsibility.
    #flux_field = get_ref_flux_field(refCat.schema, filterName=meta.filterName)

    for shard, is_on_boundary in zip(shards, boundaryMask):
        if shard is None:
            continue
        if is_on_boundary:
            refCat.extend(trim2circle(shard, circle.coord, circle.radius))
        else:
            refCat.extend(shard)
    
    # trim2bbox will not work without contiguous catalog, however deep copy
    # at this location does not seem to work for some reason 
    #if not refCat.isContiguous():
    #    refCat = refCat.copy(deep=True)

    # add and initialize centroid and hasCentroid fields after loading to avoid
    # wasting space in the saced catalogs, that is, in our loaded in-memory catalogs
    expandedCat = utils.add_centroid_columns(refCat)

    # retun Structs or break the continuity with LSST Stack?
    #pipeBase.Struct(refCat=expandedCat, fluxField=flux_field)
    return expandedCat 


def create_reference_catalogs(refCatConf, refCatLoc, refObjConf, sorted_ = True):
    """
    Wrapper function that is sets-up and returns a function that creates
    the reference catalog for an exposure from catalog shards overlapping it.

    Parameters
    ----------
    refCatConf : reference catalog configuration
    refCatLoc  : a path to the reference catalog
    refObjConf : a LoadReferenceObject configuration that defines the pairing
                 ("translation") between reference catalog and science catalog
                 schemas.
    sorted_    : returned reference catalog will be sorted by default. Useful
                 when the next step is matching to a source catalog, otherwise
                 not required.
    """
    def create_refCat(exposure):
        """
        Returns a catalog of objects in the exposure from a reference catalog.

        Parameters
        ----------
        exposure : exposureF object, or simmilar
        """

        meta = utils.get_exp_metadata(exposure)

        circle, shard_ids, boundary_mask = resolve_shard_ids(refCatConf, exposure)
        shard_locs = utils.generate_shard_paths(shard_ids, refCatConf, refCatLoc)
        shards     = utils.get_shards(shard_locs)
        
        mergeRes = merge_and_mask2circle(circle, shards, boundary_mask, 
                                         refCatConf, refCatLoc, refObjConf)

        # trim to bbox seems to break contiguity again
#        if not mergeRes.isContiguous():
#            mergeRes = mergeRes.copy(deep=True)

        refCat = trim2bbox(mergeRes, bbox=circle.bbox, wcs=meta.wcs)


        if not refCat.isSorted() and sorted_:
            refCat.sort()
        # trim to bbox and sorting break contiguity again - order important
        if not refCat.isContiguous():
            refCat = refCat.copy(deep=True)

        return refCat
    return create_refCat




################################################################################
###############################  Matchers  #####################################
################################################################################

def match_sources(wcs, onFilter):
    """
    Wrapper function that is sets-up and returns a function that matches
    sources between a catalog pair.
    
    Parameters
    ----------
    wcs      : WCS for the catalog pair in question
    onFilter : on which filter should the matching be performed
    """ 
    def src_mtch(catalogPair):
        """
        Match objects in a catalog pair.
        
        Parameters
        ----------
        catalogPair : a list or a tuple of (scienceCat, referenceCat).
        """ 
        scienceCat, referenceCat = catalogPair
        
        sciSrcSelTask = sourceSelectorRegistry["science"]()
        sciSelSrcs    = sciSrcSelTask.run(scienceCat)
        
        refSrcSelTask = sourceSelectorRegistry["references"]()
        refSelSrcs    = refSrcSelTask.run(referenceCat)
        
        matchTask = MatchPessimisticBTask()
        refFluxField = utils.get_ref_flux_field(referenceCat.getSchema(), 
                                                filterName=onFilter)
        matches = matchTask.matchObjectsToSources(refCat = refSelSrcs.sourceCat,
                                                  sourceCat = sciSelSrcs.sourceCat,
                                                  wcs = wcs,
                                                  refFluxField = refFluxField,
                                                  match_tolerance = None)

        packed_matches = afwTable.packMatches(matches.matches)            
        return packed_matches
    return src_mtch


def match_sources_with_metadata(catalogPairWithMetadata):
    """
    Match sources in a catalog pair given exposure (from which the catalogs were produced) metadata.
    
    """
    (scienceCat, referenceCat), meta = catalogPairWithMetadata
    fluxField = utils.get_ref_flux_field(referenceCat.schema, 
                                         filterName=meta.filterName)

    sciSrcSelTask = sourceSelectorRegistry["science"]()
    sciSelSrcs    = sciSrcSelTask.run(scienceCat)
    
    refSrcSelTask = sourceSelectorRegistry["references"]()
    refSelSrcs    = refSrcSelTask.run(referenceCat)
    
    matchTask = MatchPessimisticBTask()
    matches = matchTask.matchObjectsToSources(
        refCat          = refSelSrcs.sourceCat,
        sourceCat       = sciSelSrcs.sourceCat,
        wcs             = meta.wcs,
        refFluxField    = fluxField,
        match_tolerance = None,
    )

    packed_matches = afwTable.packMatches(matches.matches)
    return packed_matches
        
