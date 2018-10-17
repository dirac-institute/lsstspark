import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.geom  as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math  as afwMath
import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms as measAlgo
import lsst.daf.persistence as dafPersist

from . import utils

__all__ = ["match_and_fit_wcs", "wcs_solver"]

# matchDistanceSigma - maximum match distance is set to mean match distance
# + matchDistanceSigma*std_dev_match_distance, ignored if not fitting a WCS
def compute_match_stats_on_sky(match_list, matchDistanceSigma=2):
    distStatsInRadians = measAstrom.makeMatchStatistics(match_list, 
                                                        afwMath.MEANCLIP | afwMath.STDEVCLIP)
    # lsst.geom.radians only exists in a newer stack verion, in the old one
    # it was afw.geom.radians
    distMean = distStatsInRadians.getValue(afwMath.MEANCLIP)*afwGeom.radians
    distStdDev = distStatsInRadians.getValue(afwMath.STDEVCLIP)*afwGeom.radians
    return pipeBase.Struct(
        distMean=distMean,
        distStdDev=distStdDev,
        maxMatchDist=distMean + matchDistanceSigma*distStdDev,
    )


def match_and_fit_wcs(refCat, sourceCat, refFluxField, bbox, wcs, match_tolerance):
    matcher = measAstrom.MatchPessimisticBTask()
    matchRes = matcher.matchObjectsToSources(
        refCat = refCat, 
        sourceCat = sourceCat,
        wcs = wcs,
        refFluxField = refFluxField,
        match_tolerance = match_tolerance,
    )

    fitter = measAstrom.FitTanSipWcsTask()
    fitRes = fitter.fitWcs(
        matches = matchRes.matches,
        initWcs = wcs,
        bbox = bbox,
        refCat = refCat,
        sourceCat = sourceCat,
        exposure = None,
    )

    fitWcs = fitRes.wcs
    scatterOnSky = fitRes.scatterOnSky

    return pipeBase.Struct(
        matches=matchRes.matches,
        wcs=fitWcs,
        scatterOnSky=scatterOnSky,
        match_tolerance=matchRes.match_tolerance,
    )


def wcs_solver(max_iter=3, minMatchDistanceArcSec=0.001):
    # note that you can also have a catalogPairWithMetadata, but you can't 
    # update the exposure WCS then unless you map the returned WCS's onto 
    # ExposureRDD later on - which in itself is nto that bad, but I think
    # maybe too un-natural?
    def solver(catalogPairWithExposure):
        (sourceCat, refCat), exposure = catalogPairWithExposure
        exposureMetadata = utils.get_exp_metadata(exposure)
        fluxField = utils.get_ref_flux_field(refCat.schema, 
                                             filterName = exposureMetadata.filterName)

        match_tolerance = None
        wcs = exposureMetadata.wcs
        iterNum = 0
        for i in range(max_iter):
            iterNum += 1
            try:
                tryRes = match_and_fit_wcs(
                    refCat = refCat, # refObjLoader.loadPixelbox.refCat so freshly constructed?
                    sourceCat = sourceCat,
                    wcs = wcs,
                    refFluxField = fluxField,
                    bbox = exposureMetadata.bbox,
                    match_tolerance = match_tolerance,
                )
            except Exception as e:
                if i>0:
                    iterNum-=1
                    break
                else:
                    raise

            match_tolerance = tryRes.match_tolerance
            tryMatchDist = compute_match_stats_on_sky(tryRes.matches)
                    
            res = tryRes
            wcs = res.wcs
            maxMatchDist = tryMatchDist.maxMatchDist
            if maxMatchDist.asArcseconds() < minMatchDistanceArcSec:
                # there's a log entry here "that's good enough
                break
            match_tolerance.maxMatchDist = maxMatchDist

        # this is very confusing because self.usedKey is set to None unless
        # astrometry task is instantiated with a schema - do I lose tracking 
        # of what was used to fit WCS here?
        #for m in res.matches:
        #    if self.usedKey:
        #        m.second.set(self.usedKey, True)
        
        exposure.setWcs(res.wcs)
        return res.wcs
        
#        return pipeBase.Struct(
#            refCat = refCat,
            #matches = res.matches,
#            scatterOnSky = res.scatterOnSky,
            # I kicked this out because it's made out of refObjLoader call
            # to getMetadataBox so it should probably be with RefMatchTask?
            # I don't even know what this really is, it's like this struct
            # except all data comes straight from the exposure - can reconstruct?
            # TODO: test that there's no internal state change and that this works
#            matchMeta = matchMeta,
#        )
    return solver

        

################################################################################
#                               TESTING GROUNDS
################################################################################
#from RefMatchTask import *
#from spark_detect_sources import *
#
## get a calexp as ExposureF
#calexppath = "night_9/Night_9_Products/rerun/Night_1_rerun_ProcessCcd/0308355/calexp/calexp-0308355_01.fits"
#exposure = lsst.afw.image.ExposureF(calexppath)
#
## get metadata
#meta = utils.get_exp_metadata(exposure)
#
## create science catalog
#detect_sources = detect_and_measure(psfIterations=2, doMeasurePsf=True)
#scienceCat = detect_sources(exposure)[2]
#
## set up refObjLoader and refCat configurations
#refCatConf = DatasetConfig()
#refCatConf.ref_dataset_name  = "ps1_pv3_3pi_20170110"
#refCatConf.indexer = "HTM"
#refCatLoc = "/epyc/users/dinob/spark_proc/source_det/night_9/Night_9_Products/ref_cats"
#
## reference object setup 
#refObjConf = LoadIndexedReferenceObjectsConfig()
#refObjConf.filterMap = {'u': 'g', 'Y': 'y', "VR": "r"}
#refObjConf.ref_dataset_name = refCatConf.ref_dataset_name
#
## create referece catalog
#create_reference_catalog = create_reference_catalogs(refCatConf, refCatLoc, refObjConf)
#referenceCat = create_reference_catalog(exposure)
#
## create matches
#matches = match_sources_with_metadata([(scienceCat, referenceCat), meta])
#    
## try and run solver
#solvedRes = wcs_solver()( ((scienceCat, referenceCat), exposure) )
#print(solvedRes)

