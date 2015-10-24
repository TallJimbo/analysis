import numpy
import operator
import lsst.afw.table
from . import display

import logging
logging.basicConfig(level=logging.DEBUG)

# Prefixes for fields to take from the _ref catalog
REF_PREFIXES = ("id", "coord", "parent", "detect", "merge")

# Prefixes of fields to take from only one of _meas and _forced_src, since they're the same in both.
SHARED_PREFIXES = ("flags", "deblend", "calib")

# Flux fields that should have associated magnitude fields created for them.
MAG_FIELDS = ("flux.gaussian", "flux.psf", "flux.kron",
              "cmodel.flux", "cmodel.exp.flux", "cmodel.dev.flux", "cmodel.initial.flux")


CALCULATED_FIELDS = {}

def addRadiusFields(proxy):
    proxy._children["rDet"] = ColumnAttributeProxy(value=(proxy.xx*proxy.yy-proxy.xy*proxy.xy)**0.25)
    proxy._children["rTr"] = ColumnAttributeProxy(value=(0.5*(proxy.xx + proxy.yy))**0.5)
CALCULATED_FIELDS["ellipse"] = addRadiusFields

def addCModelFields(proxy):
    if not hasattr(proxy.exp, "ellipse"): return
    xx = (1.0-proxy.fracDev)*proxy.exp.ellipse.xx + proxy.fracDev*proxy.dev.ellipse.xx
    yy = (1.0-proxy.fracDev)*proxy.exp.ellipse.yy + proxy.fracDev*proxy.dev.ellipse.yy
    xy = (1.0-proxy.fracDev)*proxy.exp.ellipse.xy + proxy.fracDev*proxy.dev.ellipse.xy
    proxy._children["ellipse"] = ColumnAttributeProxy(
        children=dict(
            xx=ColumnAttributeProxy(value=xx),
            yy=ColumnAttributeProxy(value=yy),
            xy=ColumnAttributeProxy(value=xy),
        )
    )
    addRadiusFields(proxy._children["ellipse"])
CALCULATED_FIELDS["cmodel"] = addCModelFields

class ColumnAttributeProxy(object):

    @classmethod
    def _build(cls, columns):
        parsed = {}
        children = {}
        value = None
        for k, v in columns.iteritems():
            if len(k) == 0:
                value = v
            else:
                parsed.setdefault(k[0].replace("-", "_"), {})[k[1:]] = v
        for parsed_name, parsed_columns in parsed.iteritems():
            children[parsed_name] = ColumnAttributeProxy._build(parsed_columns)
            calculated = CALCULATED_FIELDS.get(parsed_name, None)
            if calculated:
                calculated(children[parsed_name])
        return cls(children, value)

    def __init__(self, children=None, value=None):
        if children is None: children = {}
        self.value = value
        self._children = children
        assert isinstance(self._children, dict)

    def __array__(self):
        if self.value is None:
            raise ValueError("No value associated with name")
        return self.value

    def __eq__(self, other): return operator.eq(self.value, other)
    def __ne__(self, other): return operator.ne(self.value, other)
    def __gt__(self, other): return operator.gt(self.value, other)
    def __lt__(self, other): return operator.lt(self.value, other)
    def __ge__(self, other): return operator.ge(self.value, other)
    def __le__(self, other): return operator.le(self.value, other)
    def __neg__(self): return operator.neg(self.value)
    def __pos__(self): return operator.pos(self.value)
    def __pow__(self, other): return operator.pow(self.value, other)
    def __add__(self, other): return operator.add(self.value, other)
    def __radd__(self, other): return operator.add(other, self.value)
    def __sub__(self, other): return operator.sub(self.value, other)
    def __rsub__(self, other): return operator.sub(other, self.value)
    def __mul__(self, other): return operator.mul(self.value, other)
    def __rmul__(self, other): return operator.mul(other, self.value)
    def __div__(self, other): return operator.div(self.value, other)
    def __rdiv__(self, other): return operator.div(other, self.value)
    def __truediv__(self, other): return operator.truediv(self.value, other)
    def __rtruediv__(self, other): return operator.truediv(other, self.value)
    def __floordiv__(self, other): return operator.floordiv(self.value, other)
    def __rfloordiv__(self, other): return operator.floordiv(other, self.value)

    def __dir__(self):
        names = self._children.keys()
        names.sort()
        return names

    def __getattr__(self, name):
        return self._children[name]

    def _get_str_lines(self, depth=1, key=(), **kwds):
        lines = []
        for name, child in sorted(self._children.items()):
            if depth > 0:
                lines.extend(child._get_str_lines(depth-1, key + (name,), **kwds))
            else:
                lines.append((".".join(key), "..."))
        if self.value is not None:
            s = numpy.array_str(self.value, **kwds).replace("\n", " ")
            if key:
                lines.append((".".join(key), s))
            else:
                lines.append(("value", s))
        return lines

    def __str__(self):
        return self.show()

    def show(self, depth=3, **kwds):
        lines = self._get_str_lines(depth=3, **kwds)
        max_key_length = max(len(item[0]) for item in lines) + 1
        fmt = "{{:{}s}} {{}}".format(max_key_length)
        return "\n".join(fmt.format(key + ":", value) for key, value in lines)

    def __len__(self):
        if self.value is None:
            return len(self._children.itervalues().next())
        return len(self.value)

    def __getitem__(self, k):
        return type(self)(
            children={name: child[k] for name, child in self._children.iteritems()},
            value=(self.value[k] if self.value is not None else None),
        )

    def ravel(self, *args, **kwds):
        return self.value.ravel(*args, **kwds)


class ObjectCatalog(ColumnAttributeProxy):

    @classmethod
    def read(cls, butler, dataIds=(), tracts=(), tract=None, patches=(), patch=None,
             filters=None, filter=None, forced=True, meas=True, images=True, footprints="heavy"):

        dataIds = list(dataIds)
        if filters is None:
            if filter is None:
                filters = ("g", "r", "i", "z", "y")
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
                dataIds.append(dict(tract=t, patch=p))

        refCats = []
        totalSize = 0
        for dataId in dataIds:
            logging.debug("Reading deepCoadd_ref for {}".format(dataId))
            newRefCat = butler.get("deepCoadd_ref", dataId, immediate=True,
                                   flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)
            refCats.append(newRefCat)
            totalSize += len(newRefCat)

        columns = {
            ("tract",): numpy.zeros(totalSize, dtype=int),
            ("patch",): numpy.zeros(totalSize, dtype="S5"),
            }
        if footprints:
            for b in filters:
                columns[(b, "footprint")] = numpy.zeros(totalSize, dtype=object)

        if images:
            coadds = {b: {} for b in filters}
        else:
            coadds = None

        if footprints == "heavy":
            measLoadFlags = 0
        elif footprints:
            measLoadFlags = lsst.afw.table.SOURCE_IO_NO_HEAVY_FOOTPRINTS
        else:
            measLoadFlags = lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS

        lsst.afw.image.Calib.setThrowOnNegativeFlux(False)

        offset = 0
        for n, dataId in enumerate(dataIds):
            size = len(refCats[n])

            def assignCol(key, subcol):
                if key not in columns:
                    columns[key] = numpy.zeros((totalSize,) + subcol.shape[1:], dtype=subcol.dtype)
                columns[key][offset:offset+size] = subcol

            columns[("tract",)][:] = dataId["tract"]
            columns[("patch",)][:]= dataId["patch"]

            d = refCats[n].extract("*")
            for name, subcol in d.iteritems():
                if any(name.startswith(p) for p in REF_PREFIXES):
                    assignCol(tuple(name.split(".")), subcol)
            refCats[n] = None   # allow garbage collection

            for b in filters:

                if meas or forced or images:
                    logging.debug("Reading deepCoadd for {}, {}".format(b, dataId))
                    coadd = butler.get("deepCoadd", dataId, immediate=True,
                                       filter="HSC-"+b.upper())
                    calib = coadd.getCalib()

                if meas:
                    logging.debug("Reading deepCoadd_meas for {}, {}".format(b, dataId))
                    measCat = butler.get("deepCoadd_meas", dataId, immediate=True,
                                         filter="HSC-"+b.upper(), flags=measLoadFlags)
                    d = measCat.extract("*")
                    for name, subcol in d.iteritems():
                        if any(name.startswith(p) for p in REF_PREFIXES):
                            continue
                        if any(name.startswith(p) for p in SHARED_PREFIXES):
                            key = (b,) + tuple(name.split("."))
                        else:
                            key = (b, "meas") + tuple(name.split("."))
                        assignCol(key, subcol)
                        if name in MAG_FIELDS:
                            mag, magErr = calib.getMagnitude(subcol, d[name + ".err"])
                            magKey = (b, "meas") + tuple(name.replace("flux", "mag").split("."))
                            assignCol(magKey, mag)
                            assignCol(magKey + ("err",), magErr)
                    if footprints:
                        fpCol = columns[(b, "footprint")]
                        if images:
                            logging.debug("Fixing DETECTED mask plane for {}, {}".format(b, dataId))
                            mask = coadd.getMaskedImage().getMask()
                            detPlane = mask.getMaskPlane("DETECTED")
                            detBits = mask.getPlaneBitMask("DETECTED")
                            mask.clearMaskPlane(detPlane)
                        for i, record in enumerate(measCat):
                            fpCol[i+offset] = record.getFootprint()
                            if images:
                                lsst.afw.detection.setMaskFromFootprint(mask, record.getFootprint(), detBits)
                    del measCat

                if forced:
                    logging.debug("Reading deepCoadd_meas for {}, {}".format(b, dataId))
                    forcedCat = butler.get("deepCoadd_forced_src", dataId, immediate=True,
                                           filter="HSC-" + b.upper(),
                                           flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)
                    d = forcedCat.extract("*")
                    for name, subcol in d.iteritems():
                        if any(name.startswith(p) for p in REF_PREFIXES):
                            continue
                        if any(name.startswith(p) for p in SHARED_PREFIXES):
                            if meas:
                                continue
                            else:
                                key = (b,) + tuple(name.split("."))
                        else:
                            key = (b, "forced") + tuple(name.split("."))
                        assignCol(key, subcol)
                        if name in MAG_FIELDS:
                            mag, magErr = calib.getMagnitude(subcol, d[name + ".err"])
                            magKey = (b, "forced") + tuple(name.replace("flux", "mag").split("."))
                            assignCol(magKey, mag)
                            assignCol(magKey + ("err",), magErr)
                    del forcedCat

                if images:
                    coadds[b].setdefault(dataId["tract"], {}).setdefault(dataId["patch"], coadd)

            offset += size

        logging.debug("Creating views")
        self = cls._build(columns)
        self._coadds = coadds
        self.filters = filters
        return self

    def coadd(self, filter, tract=None, patch=None):
        d1 = self._coadds[filter]
        if tract is None:
            d2 = d1.itervalues().next()
        else:
            d2 = d1[tract]
        if patch is None:
            return d2.itervalues().next()
        else:
            return d2[patch]

    def display(self, tract=None, patch=None, frames=None, frame0=0):
        return display.CoaddDisplay(self, tract, patch, frames=None, frame0=frame0)

    def __getitem__(self, k):
        r = ColumnAttributeProxy.__getitem__(self, k)
        r.coadds = self._coadds
        r.filters = self.filters
        return r


if __name__ == "__main__":
    import lsst.daf.persistence
    butler = lsst.daf.persistence.Butler("/home/jbosch/HSC/data/rerun/HSC-1339/baseline")
    objs = ObjectCatalog.read(butler, tract=8522, patch="5,5", filters=("i",), forced=False)