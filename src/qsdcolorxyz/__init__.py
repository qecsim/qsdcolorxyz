"""
**qsdcolorxyz** is a Python 3 package that extends `qecsim`_ with additional
components for the XYZ variant of the color code.

See the README at `qsdcolorxyz`_ for details.

.. _qecsim: https://github.com/qecsim/qecsim
.. _qsdcolorxyz: https://github.com/qecsim/qsdcolorxyz
"""

# import classes in dependency order
from ._color666xyzpauli import Color666XYZPauli  # noqa: F401
from ._color666xyzcode import Color666XYZCode  # noqa: F401
from ._color666xyzmpsdecoder import Color666XYZMPSDecoder  # noqa: F401

__version__ = '0.1b1'
