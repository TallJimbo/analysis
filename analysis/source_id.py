from lsst.obs.hsc import HscMapper

def splitCoaddId(oid, asDict=True, hasFilter=True):
    """Split an ObjectId (maybe an numpy array) into tract, patch, [filter], and objId.
    See obs/subaru/python/lsst/obs/hscSim/hscMapper.py
    """

    oid = np.array(oid, dtype='int64')
    objId = np.bitwise_and(oid, 2**mapper._nbit_id - 1)
    oid >>= mapper._nbit_id

    if hasFilter:
        filterId = np.bitwise_and(oid, 2**HscMapper._nbit_filter - 1).astype('int32')
        oid >>= HscMapper._nbit_filter

        filterName = np.empty(oid.size, "a6")

        if filterId.size == 1:
            filterId = [int(filterId)] # as you can't iterate over a length-1 np array

        for fid in set(filterId):
            name = afwImage.Filter(int(fid)).getName()

            filesystemName = "HSC-%s" % name.upper() # name mapper needs
            try:
                afwImage.Filter(filesystemName)
                name = filesystemName
            except:
                pass

            filterName[filterId == fid] = name
    else:
        filterName = None

    patchY = np.bitwise_and(oid, 2**HscMapper._nbit_patch - 1).astype('int32')
    oid >>= HscMapper._nbit_patch
    patchX = np.bitwise_and(oid, 2**HscMapper._nbit_patch - 1).astype('int32')
    oid >>= HscMapper._nbit_patch
    add = np.core.defchararray.add # why isn't this easier to find?
    patch = add(add(patchX.astype(str), ","), patchY.astype(str))
    patch.shape = filterName.shape # why do I have to do this?

    tract = oid.astype('int32')

    if oid.size == 1:     # sqlite doesn't like numpy types
        filterName = str(filterName[0])
        tract = int(tract)
        patch = str(patch[0])
        objId = int(objId)

    if asDict:
        return {"filter" : filterName, "tract" : tract, "patch" : patch, "objId" : objId}
    else:
        return filterName, tract, patch, objId
