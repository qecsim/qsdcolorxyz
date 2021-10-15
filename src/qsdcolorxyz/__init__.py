"""
`qsdcolorxyz`_ is a Python 3 package that extends `qecsim`_ with additional
components for the XYZ variant of the color code.

.. _qsdcolorxyz: https://github.com/qecsim/qsdcolorxyz
.. _qecsim: https://github.com/qecsim/qecsim

The key difference of the XYZ color code, relative to standard color code, is
that the single-qubit Pauli operators of the plaquette stabilizers and logical
operators are permuted X->Y->Z->X on sites identified with one color of a
bi-coloring of the color code lattice.

Compare the `plaquette`_ stabilizers and `logical`_ operators of the
`standard color code`_ with the
:meth:`permuted plaquette <qsdcolorxyz.Color666XYZPauli.plaquette>` stabilizers
and :meth:`permuted logical <qsdcolorxyz.Color666XYZPauli.logical_x>` operators
of the :class:`XYZ color code <qsdcolorxyz.Color666XYZCode>`.

.. _standard color code: https://qecsim.github.io/api/models/color.html#qecsim-models-color-color666code
.. _plaquette: https://qecsim.github.io/api/models/color.html#qecsim.models.color.Color666Pauli.plaquette
.. _logical: https://qecsim.github.io/api/models/color.html#qecsim.models.color.Color666Pauli.logical_x
"""

# import classes in dependency order
from ._color666xyzpauli import Color666XYZPauli  # noqa: F401
from ._color666xyzcode import Color666XYZCode  # noqa: F401
from ._color666xyzmpsdecoder import Color666XYZMPSDecoder  # noqa: F401

__version__ = '0.1b1'
