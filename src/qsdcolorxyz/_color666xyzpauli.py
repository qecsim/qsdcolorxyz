from qecsim.models.color import Color666Pauli


class Color666XYZPauli(Color666Pauli):
    """
    Defines a Pauli operator on a color 6.6.6 XYZ lattice.

    Notes:

    * This is a utility class used by color 6.6.6 XYZ implementations of the core models.
    * It is typically instantiated using :meth:`qsdcolorxyz.Color666XYZCode.new_pauli`
    """

    def site(self, operator, *indices, apply_permute=False):
        """
        As :meth:`Color666Pauli.site` with additional apply_permute parameter.

        :param operator: Pauli operator. One of 'I', 'X', 'Y', 'Z'.
        :type operator: str
        :param indices: Any number of indices identifying a site in the format (row, column).
        :type indices: Any number of 2-tuple of int
        :param apply_permute: Permute the operator X->Y->Z->X if the site may be permuted.
        :type apply_permute: bool
        :return: self (to allow chaining)
        :rtype: Color666Pauli
        :raises IndexError: If index is not a site index.
        """
        for index in indices:
            if apply_permute and self.code.is_site_permutable(index):
                op = {'X': 'Y', 'Y': 'Z', 'Z': 'X'}.get(operator, operator)
            else:
                op = operator
            super().site(op, index)
        return self

    def plaquette(self, operator, index):
        r"""
        Apply a plaquette operator at the given index.

        This differs from the standard Color 6.6.6 code as follows:

        'X' plaquette:
        ::

            X-Y
            |  \
            Y   X
             \  |
              X-Y

        'Y' plaquette:
        ::

            Y-Z
            |  \
            Z   Y
             \  |
              Y-Z

        'Z' plaquette:
        ::

            Z-X
            |  \
            X   Z
             \  |
              Z-X

        Notes:

        * Index is in the format (row, column).
        * Parts of plaquettes that lie outside the lattice have no effect on the lattice.

        :param operator: Pauli operator. One of 'I', 'X', 'Y', 'Z'.
        :type operator: str
        :param index: Index identifying the plaquette in the format (row, column).
        :type index: 2-tuple of int
        :return: self (to allow chaining)
        :rtype: Color666XYZPauli
        :raises IndexError: If index is not a plaquette index.
        """
        r, c = index
        # check valid index
        if not self.code.is_plaquette(index):
            raise IndexError('{} is not a plaquette index.'.format(index))
        # flip plaquette sites
        self.site(operator, (r - 1, c - 1), (r - 1, c), (r, c - 1), (r, c + 1),
                  (r + 1, c), (r + 1, c + 1), apply_permute=True)
        return self

    def logical_x(self):
        """
        Apply a logical X operator, i.e. column of alternating X, Y along leftmost sites.

        Notes:

        * The column of operators is applied to the leftmost column to allow optimisation of the MPS decoder.

        :return: self (to allow chaining)
        :rtype: Color666XYZPauli
        """
        for row in range(self.code.bound + 1):
            index = row, 0
            if self.code.is_site(index):
                self.site('X', index, apply_permute=True)
        return self

    def logical_z(self):
        """
        Apply a logical Z operator, i.e. column of alternating Z, X along leftmost sites.

        Notes:

        * The column of operators is applied to the leftmost column to allow optimisation of the MPS decoder.

        :return: self (to allow chaining)
        :rtype: Color666XYZPauli
        """
        for row in range(self.code.bound + 1):
            index = row, 0
            if self.code.is_site(index):
                self.site('Z', index, apply_permute=True)
        return self
