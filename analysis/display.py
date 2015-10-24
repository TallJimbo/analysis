import collections
import lsst.afw.geom.ellipses
import lsst.afw.display.ds9

class CoaddDisplay(object):

    def __init__(self, objs, tract=None, patch=None, frames=None, frame0=0):
        if frames is None:
            self.frames = collections.OrderedDict([(b, n + frame0) for n, b in enumerate(objs.filters)])
        elif isinstance(frames, collections.Mapping):
            self.frames = collections.OrderedDict(frames)
        elif isinstance(frames, collections.Sequence):
            self.frames = collections.OrderedDict(zip(filters, frames))
        self.objs = objs
        self.coadds = {}
        for b in objs.filters:
            self.coadds[b] = objs.coadd(b, tract=tract, patch=patch)
            self.xy0 = self.coadds[b].getXY0()

    def images(self):
        lsst.afw.display.setMaskTransparency(80.0)
        for b, frame in self.frames.iteritems():
            lsst.afw.display.ds9.mtv(self.coadds[b], frame=frame)
            lsst.afw.display.ds9Cmd("scale limits -0.5 5.0")
        lsst.afw.display.ds9Cmd("lock frame image")
        lsst.afw.display.ds9Cmd("lock scale")

    def cmodel(self, objs=None, clear=True, fitRegions=False):
        if objs is None:
            objs = self.objs
        for b, frame in self.frames.iteritems():
            if clear:
                lsst.afw.display.ds9Cmd("frame {}".format(frame))
                lsst.afw.display.ds9Cmd("regions delete all")
            eI = getattr(objs, b).meas.cmodel.initial.ellipse
            eF = getattr(objs, b).meas.cmodel.ellipse
            if fitRegions:
                eRI = getattr(objs, b).meas.cmodel.region.initial.ellipse
                eRF = getattr(objs, b).meas.cmodel.region.final.ellipse
            x = getattr(objs, b).meas.centroid.sdss.x.value
            y = getattr(objs, b).meas.centroid.sdss.y.value
            with lsst.afw.display.Buffering():
                for i in xrange(len(objs)):
                    qI = lsst.afw.geom.ellipses.Quadrupole(eI.xx.value[i], eI.yy.value[i], eI.xy.value[i])
                    qF = lsst.afw.geom.ellipses.Quadrupole(eF.xx.value[i], eF.yy.value[i], eF.xy.value[i])
                    lsst.afw.display.ds9.dot(
                        qI, x[i]-self.xy0.getX(), y[i]-self.xy0.getY(),
                        ctype='yellow', frame=frame
                    )
                    lsst.afw.display.ds9.dot(
                        qF, x[i]-self.xy0.getX(), y[i]-self.xy0.getY(),
                        ctype='green', frame=frame
                    )
                    if fitRegions:
                        qRI = lsst.afw.geom.ellipses.Quadrupole(eRI.xx.value[i], eRI.yy.value[i], eRI.xy.value[i])
                        qRF = lsst.afw.geom.ellipses.Quadrupole(eRF.xx.value[i], eRF.yy.value[i], eRF.xy.value[i])
                        lsst.afw.display.ds9.dot(
                            qRI, x[i]-self.xy0.getX(), y[i]-self.xy0.getY(),
                            ctype='magenta', frame=frame
                        )
                        lsst.afw.display.ds9.dot(
                            qRF, x[i]-self.xy0.getX(), y[i]-self.xy0.getY(),
                            ctype='blue', frame=frame
                        )


    def footprints(self, objs=None, clear=True):
        if objs is None:
            objs = self.objs
        for b, frame in self.frames.iteritems():
            if clear:
                lsst.afw.display.ds9Cmd("frame {}".format(frame))
                lsst.afw.display.ds9Cmd("regions delete all")
            fps = getattr(objs, b).footprint
            with lsst.afw.display.Buffering():
                for i in xrange(len(objs)):
                    lsst.afw.display.utils.drawFootprint(fps.value[i], XY0=self.xy0, ctype="cyan",
                                                         peaks=True, ctypePeak="cyan", frame=frame)


    def kron(self, objs=None, clear=True, radiusForRadius=False):
        if objs is None:
            objs = self.objs
        for b, frame in self.frames.iteritems():
            if clear:
                lsst.afw.display.ds9Cmd("frame {}".format(frame))
                lsst.afw.display.ds9Cmd("regions delete all")
            shape = getattr(objs, b).meas.shape.sdss
            r1 = getattr(objs, b).meas.flux.kron.radius.value
            r2 = getattr(objs, b).meas.flux.kron.radiusForRadius.value
            x = getattr(objs, b).meas.centroid.sdss.x.value
            y = getattr(objs, b).meas.centroid.sdss.y.value
            with lsst.afw.display.Buffering():
                for i in xrange(len(objs)):
                    ellipse = lsst.afw.geom.ellipses.Quadrupole(
                        shape.xx.value[i],
                        shape.yy.value[i],
                        shape.xy.value[i]
                    )
                    r0 = ellipse.getDeterminantRadius()
                    ellipse.scale(float(r1[i]/r0))
                    lsst.afw.display.ds9.dot(
                        ellipse, x[i]-self.xy0.getX(), y[i]-self.xy0.getY(),
                        ctype='orange', frame=frame
                    )
                    if radiusForRadius:
                        ellipse.scale(float(r2[i]/r1[i]))
                        lsst.afw.display.ds9.dot(
                            ellipse, x[i]-self.xy0.getX(), y[i]-self.xy0.getY(),
                            ctype='red', frame=frame
                        )

    def annotate(self, objs):
        if objs is None:
            objs = self.objs
        for b, frame in self.frames.iteritems():
            psfMag = getattr(objs, b).meas.mag.psf.value
            cmodelMag = getattr(objs, b).meas.cmodel.mag.value
            x = getattr(objs, b).meas.centroid.sdss.x.value
            y = getattr(objs, b).meas.centroid.sdss.y.value
            for i in xrange(len(objs)):
                with lsst.afw.display.Buffering():
                    lsst.afw.display.ds9.dot(
                        "    psf={psf}, cmodel={cmodel}".format(psf=psfMag[i], cmodel=cmodelMag[i]),
                        x[i] - self.xy0.getX(), y[i] - self.xy0.getY(),
                        ctype='red', frame=frame
                    )