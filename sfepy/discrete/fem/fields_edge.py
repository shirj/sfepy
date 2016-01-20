"""
"""
import numpy as nm

from sfepy.base.base import assert_, basestr, Struct
from sfepy.discrete.common.fields import parse_shape
from sfepy.discrete.fem.fields_base import VolumeField
from sfepy.discrete.fem.poly_spaces import PolySpace


class HDivField(VolumeField):
    family_name = 'volume_Hdiv_Raviart-Thomas_discontinuous'

    def __init__(self, name, dtype, shape, region, approx_order=1):
        """
        Create a finite element field.

        Parameters
        ----------
        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or (2,)
            or 3 or (3,)) or 'vector', or a tuple. The field shape determines
            the shape of the FE base functions and is related to the number of
            components of variables and to the DOF per node count, depending
            on the field kind.
        region : Region
            The region where the field is defined.
        approx_order : int or tuple
            The FE approximation order. The tuple form is (order, has_bubble),
            e.g. (1, True) means order 1 with a bubble function.

        Notes
        -----
        Assumes one cell type for the whole region!
        """
        shape = parse_shape(shape, region.domain.shape.dim)
        if not self._check_region(region):
            raise ValueError('unsuitable region for field %s! (%s)' %
                             (name, region.name))
        assert_((len(shape) == 1) and (shape[0] == region.dim))

        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)
        self.domain = self.region.domain
        self._set_approx_order(approx_order)
        self._setup_geometry()
        self._setup_kind()
        self._setup_shape()

        self._create_interpolant()
        from debug import debug; debug()
        self._setup_global_base()
        self.setup_coors()
        self.clear_mappings(clear_all=True)

    def _set_approx_order(self, approx_order):
        """
        Set a uniform approximation order.
        """
        if isinstance(approx_order, tuple):
            self.approx_order = approx_order[0]
            self.force_bubble = approx_order[1]

        else:
            self.approx_order = approx_order
            self.force_bubble = False

        if (self.approx_order != 0) or self.force_bubble:
            raise NotImplementedError('only RT0 basis is supported!')

    def _setup_shape(self):
        """
        Setup the field's shape-related attributes, see :class:`Field`.
        """
        self.n_components = 1
        self.val_shape = self.shape

    def _create_interpolant(self):
        name = '%s_%s_%s_%d' % (self.gel.name, self.space,
                                  self.poly_space_base, self.approx_order)
        self.poly_space = PolySpace.any_from_args(name, self.gel,
                                                  self.approx_order,
                                                  base=self.poly_space_base)

        cc = nm.array([[0, 0], [0.1, 0.1]])
        print self.poly_space.eval_base(cc)
        print self.poly_space.eval_base(cc, diff=1)
        from debug import debug; debug()

    def _setup_global_base(self):
        """
        Setup global DOF/base functions, their indices and connectivity of the
        field.
        """
        n_cell = self.region.get_n_cells(True)

        self.econn = nm.zeros((n_cell, self.poly_space.n_nod), dtype=nm.int32)
        self.ori = nm.zeros_like(self.econn)
