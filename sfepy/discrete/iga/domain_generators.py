"""
IGA domain generators.
"""
import numpy as nm

from sfepy.base.base import output, Struct
import sfepy.discrete.iga as iga
from sfepy.discrete.iga.domain import NurbsPatch

def gen_patch_block_domain(dims, shape, centre, degrees, continuity=None,
                           name='block', verbose=True):
    """
    Generate a single IGA patch block in 2D or 3D of given degrees and
    continuity using igakit.

    Parameters
    ----------
    dims : array of D floats
        Dimensions of the block.
    shape : array of D ints
        Numbers of unique knot values along each axis.
    centre : array of D floats
        Centre of the block.
    degrees : array of D floats
        NURBS degrees along each axis.
    continuity : array of D ints, optional
        NURBS continuity along each axis. If None, `degrees-1` is used.
    name : string
        Domain name.
    verbose : bool
        If True, report progress of the domain generation.

    Returns
    -------
    nurbs : NurbsPatch instance
        The NURBS data. The igakit NURBS object is stored as `nurbs` attribute.
    bmesh : Struct instance
        The Bezier mesh data.
    regions : dict
        The patch surface regions.
    """
    import igakit.cad as cad

    dims = nm.asarray(dims, dtype=nm.float64)
    shape = nm.asarray(shape, dtype=nm.int32)
    centre = nm.asarray(centre, dtype=nm.float64)
    degrees = nm.asarray(degrees, dtype=nm.int32)

    if continuity is None:
        continuity = degrees - 1

    else:
        continuity = nm.asarray(continuity, dtype=nm.int32)

    dim = len(shape)

    output('generating NURBS...', verbose=verbose)

    dd = centre - 0.5 * dims
    block = cad.grid(shape - 1, degree=degrees, continuity=continuity)

    for ia in xrange(dim):
        block.scale(dims[ia], ia)

    for ia in xrange(dim):
        block.translate(dd[ia], ia)

    # Force uniform control points. This prevents coarser resolution inside the
    # block.
    shape = nm.asarray(block.points.shape[:-1])
    n_nod = nm.prod(shape)
    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    ngrid = nm.mgrid[[slice(ii) for ii in shape]]
    ngrid.shape = (dim, n_nod)

    coors = x0 + ngrid.T * dd
    coors.shape = shape.tolist() + [dim]

    ####
    blow = 1.0
    delta = 0.25 * dims[0]
    power = 2.0

    cc = coors.reshape((-1, dim))

    dirs = (cc - centre[None, :])

    radii = nm.abs(dirs).max(1)
    ii = nm.abs(radii > 1e-12)
    dirs[ii] = dirs[ii] / radii[ii][:, None]

    dirs2 = (cc - centre[None, :])
    radii2 = nm.sqrt((dirs2 * dirs2).sum(1))
    ii = nm.abs(radii2 > 1e-12)
    dirs2[ii] = dirs2[ii] / radii2[ii][:, None]

    mod = lambda x: nm.where(x >= delta, delta + ((nm.maximum(x, delta) - delta) / (x.max() - delta))**power, x)

    ## cc = cc + blow * mod[:, None] * dirs2

    cc = centre[None, :] + (blow
                            * (dirs2 * mod(radii)[:, None] * nm.sqrt(2)))

    ## cc = centre[None, :] + ((dirs2 * radii[:, None] * nm.sqrt(2)))

    ## fac = dims[0] / (cc[:, 0].max() - cc[:, 0].min())
    ## cc *= fac

    coors = cc.reshape(coors.shape)
    ####

    block.array[..., :dim] = coors

    output('...done', verbose=verbose)

    # Compute Bezier extraction data.
    output('computing Bezier mesh...', verbose=verbose)
    cs = iga.compute_bezier_extraction(block.knots, block.degree)
    n_els = [len(ii) for ii in cs]
    conn, bconn = iga.create_connectivity(n_els, block.knots, block.degree)

    ccs = iga.combine_bezier_extraction(cs)

    cps = block.points[..., :dim].copy()
    cps = cps.reshape((-1, dim))
    bcps, bweights = iga.compute_bezier_control(cps, block.weights.ravel(), ccs,
                                                conn, bconn)

    nurbs = NurbsPatch(block.knots, degrees, cps, block.weights.ravel(), cs,
                       conn)
    nurbs.nurbs = block
    bmesh = Struct(name='bmesh', cps=bcps, weights=bweights, conn=bconn)
    output('...done', verbose=verbose)

    output('defining regions...', verbose=verbose)
    regions = iga.get_patch_box_regions(n_els, block.degree)
    output('...done', verbose=verbose)

    return nurbs, bmesh, regions
