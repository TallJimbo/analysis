#!/usr/bin/env python

import uuid
import itertools
import os
import re

import numpy
import matplotlib.colors
import matplotlib.cm

import lsst.afw.geom.ellipses
import lsst.afw.display

def rgba2hex(rgba):
    return "#{0:02x}{1:02x}{2:02x}".format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

def property_setter(func, name=None):
    if name is None:
        name = func.__name__
    def getter(self):
        return self._properties[name]
    def setter(self, value):
        self._properties[name] = func(self, value)
    return property(getter, setter)

def make_extractor_property(name, allow_none=False, allow_str=False):
    def setter(self, value):
        if value is None:
            if allow_none:
                return None
            else:
                raise ValueError("Property '{}' cannot be None".format(name))
        if isinstance(value, basestring):
            key = self.catalog.schema.find(value).key
            return lambda r: r.get(key)
        else:
            return value
    return property_setter(setter, name)

def make_bool_property(name):
    def getter(self):
        return self._properties[name]
    def setter(self, value):
        assert bool(value) == value
        self._properties[name] = bool(value)


class CatalogView(object):


    def __init__(self, catalog, frame=0, tag=None, coordsys="wcsa", **kwds):
        self.catalog = catalog
        self.frame = frame
        self.tag = tag if tag is not None else uuid.uuid4().hex
        self.coordsys = coordsys
        args = dict(
            centroid=self.catalog.getCentroidDefinition(),
            ellipse=None, color="green", dash=False, text=None, width=1,
            symbol="cross", id="id", can_edit=False, can_move=False, can_rotate=False, can_delete=True
        )
        args.update(kwds)
        self._properties = {}
        self.configure(**args)

    centroid = make_extractor_property("centroid")
    id = make_extractor_property("id")
    text = make_extractor_property("text", allow_none=True)
    ellipse = make_extractor_property("ellipse", allow_none=True)
    dash = make_bool_property("dash")
    can_edit = make_bool_property("can_edit")
    can_move = make_bool_property("can_move")
    can_rotate = make_bool_property("can_rotate")
    can_delete = make_bool_property("can_delete")

    @property_setter
    def width(self, value):
        assert int(value) == value and value > 0 and value < 5
        return int(value)

    @property_setter
    def symbol(self, value):
        assert value in (None, "cross", "circle", "box", "diamond", "x")
        return value

    def configure(self, color=None, cmap=None, norm=None, vmin=None, vmax=None, **kwds):
        for name, value in kwds.iteritems():
            setattr(self, name, value)
        if color is not None:
            try:
                key = self.catalog.schema.find(color).key
                color = lambda catalog: catalog.get(key)
            except:
                if isinstance(color, basestring):
                    self._color = color
                    return
            if norm is None:
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            if cmap is None:
                cmap = matplotlib.cm.get_cmap()
            self._color = lambda catalog: cmap(norm(color(catalog)))

    def show(self):
        fmt1 = ("{s.coordsys}; "
                    "point {x} {y} # point={s.symbol} width={s.width:d} dash={s.dash:1d} "
                    "edit={s.can_edit:d} move={s.can_move:d} "
                    "rotate={s.can_rotate:d} delete={s.can_delete:d} "
                    "tag={s.tag} {text} color={{{color}}} id={id}"
                ";\n")
        fmt2 = ("{s.coordsys}; "
                    "ellipse({x},{y},{a},{b},{theta}) # width={s.width:d} dash={s.dash:1d} "
                    "edit={s.can_edit:d} move={s.can_move:d} "
                    "rotate={s.can_rotate:d} delete={s.can_delete:d} "
                    "tag={s.tag} {text} color={{{color}}} id={id}"
                ";\n")
        if isinstance(self._color, basestring):
            colors = itertools.repeat(self._color)
        else:
            colors = itertools.imap(rgba2hex, self._color(self.catalog))
        lsst.afw.display.ds9.ds9Cmd(lsst.afw.display.ds9.selectFrame(self.frame))
        xpa_cmd = "xpaset {0} regions".format(
            lsst.afw.display.ds9.getXpaAccessPoint(), self.coordsys
        )
        pfd = os.popen(xpa_cmd, "w")
        try:
            for record, color in itertools.izip(self.catalog, colors):
                x, y = self.centroid(record)
                if not numpy.isfinite(x) or not numpy.isfinite(y):
                    continue
                id = self.id(record)
                text = 'text="{}"'.format(self.text(record)) if self.text is not None else ""
                if self.symbol is not None:
                    line = fmt1.format(s=self, x=x, y=y, color=color, text=text, id=id)
                    pfd.write(line)
                    text = ""   # if we're plotting a symbol and an ellipse, only label one of them
                if self.ellipse is not None:
                    axes = lsst.afw.geom.ellipses.Axes(self.ellipse(record))
                    a = axes.getA()
                    b = axes.getB()
                    theta = lsst.afw.geom.radToDeg(axes.getTheta())
                    if not numpy.isfinite(a) or not numpy.isfinite(b) or not numpy.isfinite(theta):
                        continue
                    line = fmt2.format(s=self, x=x, y=y, color=color, text=text, id=id,
                                       a=axes.getA(), b=axes.getB(), theta=theta)
                    pfd.write(line)
        finally:
            pfd.close()

    def hide(self):
        lsst.afw.display.ds9Cmd("regions group {} delete".format(self.tag))

    def select(self):
        lsst.afw.display.ds9Cmd("regions group {} select".format(self.tag))

    def inspect(self):
        regex = re.compile("id=(\d+)")
        records = []
        for line in lsst.afw.display.ds9Cmd("regions selected".format(self.tag), get=True).split("\n"):
            m = regex.search(line)
            if m:
                print "yes: ", line
                records.append(self.catalog.find(int(m.group(1))))
            else:
                print "no: ", line
        if len(records) == 1:
            return records[0]
        elif len(records) == 0:
            return None
        else:
            catalog = type(self.catalog)(self.catalog.table)
            catalog.extend(records)
            return catalog

    def __del__(self):
        self.hide()


if __name__ == "__main__":
    import lsst.afw.table
    import lsst.afw.image
    catalog = lsst.afw.table.SourceCatalog.readFits("/home/jbosch/HSC/data/meas-HSC-I-8766-4,4.fits")
    exposure = lsst.afw.image.ExposureF("/home/jbosch/HSC/data/4,4.fits")
    #lsst.afw.display.ds9.mtv(exposure, frame=1)
    calib = exposure.getCalib()
    calib.setThrowOnNegativeFlux(False)
    mag = lambda c: calib.getMagnitude(c.get("flux.psf"))
    mask = catalog.get("flux.psf") > 0
    filtered = catalog[mask].copy(deep=True)
    cv = CatalogView(filtered, frame=1, symbol="x", ellipse="shape.sdss")
    cv.configure(color=mag, vmin=18, vmax=26)
    cv.show()
    raw_input()
    print cv.inspect().extract("*")
    cv.hide()