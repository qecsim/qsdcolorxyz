from qecsim.model import cli_description
from qecsim.models.color import Color666Code

from . import Color666XYZPauli


@cli_description('Color 6.6.6 XYZ (size INT odd >=3)')
class Color666XYZCode(Color666Code):
    r"""
    Implements a triangular color 6.6.6 XYZ code.
    """

    @property
    def label(self):
        """See :meth:`qecsim.model.StabilizerCode.label`"""
        return 'Color 6.6.6 XYZ {}'.format(self.size)

    def is_site_permutable(self, index):
        """
        Return true if the site may be permuted X->Y->Z->X.

        :param index: Index identifying a site in the format (row, column).
        :type index: 2-tuple of int
        :return: If the site may be permuted.
        :rtype: bool
        :raises IndexError: If index is not a plaquette index.
        """
        if not self.is_site(index):
            raise IndexError('{} is not a site index.'.format(index))
        return sum(index) % 3 == 1

    def new_pauli(self, bsf=None):
        """
        Convenience constructor of color 6.6.6 XYZ Pauli for this code.

        Notes:

        * For performance reasons, the new Pauli is a view of the given bsf. Modifying one will modify the other.

        :param bsf: Binary symplectic representation of Pauli. (Optional. Defaults to identity.)
        :type bsf: numpy.array (1d)
        :return: Color 6.6.6 XYZ Pauli
        :rtype: Color666XYZPauli
        """
        return Color666XYZPauli(self, bsf)
