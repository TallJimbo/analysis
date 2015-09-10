import numpy
import lsst.afw.table


class MagConverter(object):

    def __init__(self, flux, schemaMapper, mag=None, prefix=None):
        if mag is None:
            mag = flux.replace("flux", "mag")
        if prefix is not None:
            mag = "%s.%s" % (prefix, mag)
        self.fluxKey = schemaMapper.getInputSchema().find(flux).key
        self.fluxErrKey = schemaMapper.getInputSchema().find(flux + ".err").key
        self.magKey = schemaMapper.editOutputSchema().addField(mag, type=float, doc="magnitude for %s" % flux)
        self.magErrKey = schemaMapper.editOutputSchema().addField(mag + ".err", type=float,
                                                                  doc="magnitude err for %s" % flux)

    def __call__(self, catalog, calib):
        oldThrow = calib.getThrowOnNegativeFlux()
        calib.setThrowOnNegativeFlux(False)
        mag, magErr = calib.getMagnitude(catalog[self.fluxKey], catalog[self.fluxErrKey])
        catalog[self.magKey][:] = mag
        catalog[self.magErrKey][:] = magErr
        calib.setThrowOnNegativeFlux(oldThrow)


class CatalogLoader(object):

    # Prefixes for fields to take from the _ref catalog
    REF_PREFIXES = ("detect", "merge")

    # Prefixes of fields to take from only one of _meas and _forced_src, since they're the same in both.
    SHARED_PREFIXES = ("flags", "deblend", "calib")

    # Flux fields that should have associated magnitude fields created for them.
    MAG_FIELDS = ("flux.gaussian", "flux.psf", "flux.kron",
                  "cmodel.flux", "cmodel.exp.flux", "cmodel.dev.flux")

    def __init__(self, butler, filters=("g", "r", "i", "z", "y"), forced=True, meas=True):
        self.butler = butler
        self.filters = tuple(filters)
        refSchema = butler.get("deepCoadd_ref_schema", immediate=True)
        self.refMapper = lsst.afw.table.SchemaMapper(refSchema)
        self.refMapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema(), True)
        sharedFieldsAlreadyMapped = set(lsst.afw.table.SourceTable.makeMinimalSchema().getNames())
        for item in refSchema:
            prefix = item.field.getName().split(".")[0]
            if prefix in self.REF_PREFIXES:
                self.refMapper.addMapping(item.key)
                sharedFieldsAlreadyMapped.add(item.field.getName())
        self.outSchema = self.refMapper.getOutputSchema()
        self.tractKey = self.outSchema.addField("tract", type=int, doc="Coadd tract")
        self.patchXKey = self.outSchema.addField("patch.x", type=int, doc="Coadd patch X")
        self.patchYKey = self.outSchema.addField("patch.y", type=int, doc="Coadd patch Y")
        self.measMappers = dict()
        self.measMags = dict()
        if meas:
            measSchema = butler.get("deepCoadd_meas_schema", immediate=True)
            for b in self.filters:
                        continue
                self.measMappers[b] = lsst.afw.table.SchemaMapper(measSchema, self.outSchema)
                for item in measSchema:
                    if item.field.getName() in sharedFieldsAlreadyMapped:
                    prefix = item.field.getName().split(".")[0]
                    if prefix in self.SHARED_PREFIXES:
                        self.measMappers[b].addMapping(item.key, "%s.%s" % (b, item.field.getName()))
                        sharedFieldsAlreadyMapped.add(item.field.getName())
                    else:
                        self.measMappers[b].addMapping(item.key, "meas.%s.%s" % (b, item.field.getName()))
                self.measMags[b] = [MagConverter(flux, self.measMappers[b], prefix="meas.%s" % b)
                                    for flux in self.MAG_FIELDS]
                self.outSchema = measMappers[b].getOutputSchema()
        self.forcedMappers = dict()
        self.forcedMags = dict()
        if forced:
            forcedSchema = butler.get("deepCoadd_forced_src_schema", immediate=True)
            for b in self.filters:
                self.forcedMappers[b] = lsst.afw.table.SchemaMapper(forcedSchema, self.outSchema)
                for item in forcedSchema:
                    if item.field.getName() in sharedFieldsAlreadyMapped:
                        continue
                    prefix = item.field.getName().split(".")[0]
                    if prefix in self.SHARED_PREFIXES:
                        self.forcedMappers[b].addMapping(item.key, "%s.%s" % (b, item.field.getName()))
                    else:
                        self.forcedMappers[b].addMapping(item.key, "forced.%s.%s" % (b, item.field.getName()))
                self.forcedMags[b] = [MagConverter(flux, self.forcedMappers[b], prefix="forced.%s" % b)
                                      for flux in self.MAG_FIELDS]
                self.outSchema = forcedMappers[b].getOutputSchema()

    def read(self, *args, tracts=(), tract=None, patches=(), patch=None, filters=None, filter=None,
             extend=None, copy=True):
        dataIDs = list(args)
        if filters is None:
            if filter is None:
                filters = self.filters
            else:
                filters = (filter,)
        if tract is not None:
            tracts = tuple(tracts) + (tract,)
        if patch is not None:
            patches = tuple(patches) + (patch,)
        for t in tracts:
            for p in patches:
                if not isinstance(p, basestring):
                    p = "%s,%s" % p
                dataIDs.append(dict(tract=t, patch=p))
        if extend is not None:
            assert extend.schema == self.outSchema
            catalog = lsst.afw.table.SourceCatalog(extend.table)
            catalog.extend(extend)
        else:
            catalog = lsst.afw.table.SourceCatalog(self.outSchema)
        for dataID in dataIDs:
            subCat = lsst.afw.table.SourceCatalog(catalog.table)
            refCat = self.butler.get("deepCoadd_ref", immediate=True,
                                     flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS, **dataID)
            subCat.extend(refCat, self.refMapper)
            subCat[self.tractKey][:] = dataID["tract"]
            patchX, patchY = int(p) for p in dataID["patch"].split(",")
            subCat[self.patchXKey][:] = patchX
            subCat[self.patchXKey][:] = patchY
            for b, mapper in self.measMappers.iteritems():
                measCat = self.butler.get("deepCoadd_meas", immediate=True,
                                          flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS,
                                          filter="HSC-%s" % b.upper(), **dataID)
                for inRecord, outRecord in zip(measCat, subCat):
                    outRecord.assign(inRecord, mapper)
            for b, mapper in self.forcedMappers.iteritems():
                forcedCat = self.butler.get("deepCoadd_forced_src", immediate=True,
                                          flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS,
                                          filter="HSC-%s" % b.upper(), **dataID)
                for inRecord, outRecord in zip(forcedCat, subCat):
                    outRecord.assign(inRecord, mapper)
            catalog.extend(subCat)
        if copy:
            catalog = catalog.copy(deep=True)
        return catalog
