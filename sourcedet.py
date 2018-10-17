from lsst.pipe.tasks.measurePsf import MeasurePsfTask, MeasurePsfConfig
from lsst.pipe.tasks.characterizeImage import CharacterizeImageConfig
from lsst.pipe.tasks.repair import RepairTask, RepairConfig
from lsst.meas.algorithms.subtractBackground import SubtractBackgroundConfig, SubtractBackgroundTask
from lsst.meas.algorithms.detection import SourceDetectionTask, SourceDetectionConfig
from lsst.meas.base import SingleFrameMeasurementTask, SingleFrameMeasurementConfig
from lsst.afw.table import IdFactory, SourceTable
import lsst.daf.base as dafBase
from lsst.obs.base import ExposureIdInfo
import lsst.pipe.base as pipeBase

__all__ = ["detect_and_measure"]


def measure_psf(exposure, exposureIdInfo, background):
    """Performs a single estimate of the PSF for an exposure.

    Parameters
    ----------
    exposure   : an exposure of the sky (lsst.afw.ExposureF object)
    background : estimated background of exposure, a list-like object containing 
                 tuples of (afwMath.Background, interpStyle, undersampleStyle)
    exposureIdInfo : Exposure ID and number of bits used (lsst.obs.base.ExposureIdInfo)

    Returns
    ----------
    exposure : original exposure on which PSF was measured
    sources  : a catalog of detected sources (lsst.afw.table.SourceCatalog object)
               used to estimate the PSF, contains flags for objects which were used
               to estimate the PSF
    cellSet  : the spatial cell set used to determine the PSF 
               (lsst.afw.math.SpatialCellSet)
    """

    # make an empty schema - this is used to create a table. It is very
    # important that the schema is completely instantiated before
    # running any Tasks
    schema      = SourceTable.makeMinimalSchema()
    algMetadata = dafBase.PropertyList()

    # Declare all the tasks that will be run on the image - this must(!)
    # be done prior any Task is run because they silently expand the schema
    repairConf = RepairConfig()
    repairTask = RepairTask()
    
    singleFrameMeasConf = SingleFrameMeasurementConfig()
    singleFrameMeasTask = SingleFrameMeasurementTask(schema=schema, 
                                                     algMetadata=algMetadata)

    detectConf = SourceDetectionConfig()
    detectTask = SourceDetectionTask()

    measurePsfConf = MeasurePsfConfig()
    measurePsfTask = MeasurePsfTask(schema=schema) 

    # now that schema has been altered by task inits we can create a 
    # table that will be filled in with data once Tasks are run
    srcIdFactory = IdFactory.makeSource(exposureIdInfo.expId,
                                        exposureIdInfo.unusedBits)
    table        = SourceTable.make(schema, srcIdFactory)
    table.setMetadata(algMetadata)

    # It is very important Tasks are run in this very specific order
    repairTask.run(exposure, keepCRs=True)
    detRes = detectTask.run(exposure=exposure, table=table, doSmooth=True)
    singleFrameMeasTask.run(measCat=detRes.sources, exposure=exposure,
                            exposureId=exposureIdInfo.expId)
    measPsfRes = measurePsfTask.run(exposure=exposure, sources=detRes.sources, 
                                    matches=None, expId=exposureIdInfo.expId)

    # return as struct, or not?
    res = pipeBase.Struct(exposure=exposure, sources=detRes.sources, 
                          psfcells=measPsfRes.cellSet)
    return res


def detect_and_measure(psfIterations=2, doMeasurePsf=True):
    """A wrapper that returns a function that can be mapped to RDDs. 

    Parameters
    ----------
    psfIterations : number of times PSF is estimated on the exposure.
    doMeasurePsf  : boolean value, if set to False PSF will be estimated 
                    1 time. ****TO DO:****
                                      implement installing of a default Gaussian
                                      PSF
    Returns
    ----------
    Function that estimates background, PSF, repairs exposure defects, repairs CRs, 
    detects and measures sources on a given exposure.
    """
    def det_meas(exposure):
        """Estimate background and PSF, repair exposure defects and CRs, detect and
           measure sources. 

        Parameters
        ----------
        exposurePath : path to an exposure of the sky (lsst.afw.ExposureF object)

        Returns
        ----------
        exposure : exposure of the sky given by exposurePath (lsst.afw.ExposureF)
        calexp   : exposure of the sky with subtracted background, repaired defects
                   and CRs (lsst.afw.ExposureF)
        sources  : catalog of detected sources on the calexp
                   (lsst.afw.table.SourceCatalog)
        backg    : estimated background of exposure, a list-like object containing 
                   tuples of (afwMath.Background, interpStyle, undersampleStyle)
        """
        exposureIdInfo = ExposureIdInfo()
        
        bckgConfig = SubtractBackgroundConfig()
        bckgTask   = SubtractBackgroundTask()
        bckg       = bckgTask.run(exposure)
        
        NPsfIter = psfIterations if doMeasurePsf else 1
        for i in range(NPsfIter):
            dmeRes = measure_psf(exposure, exposureIdInfo=exposureIdInfo, 
                                 background=bckg)
            
        # Notice a different exp is being repaired than in the DME PSF step
        repairConf = RepairConfig()
        repairTask = RepairTask()
        repairTask.run(exposure=dmeRes.exposure)

        schema = SourceTable.makeMinimalSchema()        
        singleFrameMeasConf = SingleFrameMeasurementConfig()
        singleFrameMeasTask = SingleFrameMeasurementTask(schema=schema)
        singleFrameMeasTask.run(measCat=dmeRes.sources, exposure=dmeRes.exposure, 
                                exposureId=exposureIdInfo.expId)

        # schema can't be pickled! PSF PCA object can't be pickled!
        return (exposure, dmeRes.exposure, dmeRes.sources, bckg)
    return det_meas




#print("start!")
#calexppath = "night_9/Night_9_Products/rerun/Night_1_rerun_ProcessCcd/0308355/calexp/calexp-0308355_01.fits"
#sciexp, calexp, srcs, bckgs = detect_and_measure(psfIterations=1)(calexppath)
#print(srcs.schema)
