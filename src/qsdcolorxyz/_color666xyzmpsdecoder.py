import functools
import logging

import numpy as np
from qecsim import paulitools as pt
from qecsim import tensortools as tt
from qecsim.model import cli_description
from qecsim.models.color import Color666MPSDecoder

logger = logging.getLogger(__name__)


@cli_description('MPS ([chi] INT, ...)')
class Color666XYZMPSDecoder(Color666MPSDecoder):
    r"""
    Implements a color 6.6.6 XYZ Matrix Product State (MPS) decoder.
    """

    @classmethod
    def sample_recovery(cls, code, syndrome):
        """
        Return a sample Pauli consistent with the syndrome, created by applying a path between each plaquette identified
        by the syndrome and the nearest boundary of the same type as the plaquette.

        :param code: Color 666 XYZ code.
        :type code: Color666XYZCode
        :param syndrome: Syndrome as binary vector.
        :type syndrome: numpy.array (1d)
        :return: Sample recovery operation as color 666 XYZ Pauli.
        :rtype: Color666XYZPauli
        """
        # prepare sample
        sample_recovery = code.new_pauli()
        # iterate syndrome_indices and corrective operator
        for syndrome_indices, op in zip(code.syndrome_to_plaquette_indices(syndrome), ('Z', 'X')):
            # add path of op from syndrome_index:(r1, c1) to virtual_index:(r2, c2) inclusive
            for r1, c1 in syndrome_indices:
                # find nearest off-boundary plaquette
                r2, c2 = code.virtual_plaquette_index((r1, c1))
                # path along horizontal
                step = np.sign(c2 - c1)
                if step:
                    for cc in range(c1, c2 + step, step):
                        if code.is_site((r1, cc)):
                            sample_recovery.site(op, (r1, cc), apply_permute=True)
                # path along vertical
                step = np.sign(r2 - r1)
                if step:
                    for rr in range(r1, r2 + step, step):
                        if code.is_site((rr, c1)):
                            sample_recovery.site(op, (rr, c1), apply_permute=True)
        # return sample
        return sample_recovery

    @property
    def label(self):
        """See :meth:`qecsim.model.Decoder.label`"""
        params = [('chi', self._chi), ('tol', self._tol), ]
        return 'Color 6.6.6 XYZ MPS ({})'.format(', '.join('{}={}'.format(k, v) for k, v in params if v))

    class TNC(Color666MPSDecoder.TNC):
        """Tensor network creator"""

        @functools.lru_cache()
        def q_node_value(self, prob_dist, f, n, e, s, w, permute):
            """Return qubit tensor element value."""
            paulis = ('I', 'X', 'Y', 'Z')
            op_to_pr = dict(zip(paulis, prob_dist))
            f = pt.pauli_to_bsf(f)
            # n, e, s, w are in {0, 1, 2, 3} so create dict from index to op
            stabilizer_paulis = paulis if not permute else ('I', 'Y', 'Z', 'X')
            index_to_op = dict(zip((0, 1, 2, 3), pt.pauli_to_bsf(stabilizer_paulis)))
            # apply ops from indices to f
            op = (f + index_to_op[n] + index_to_op[e] + index_to_op[s] + index_to_op[w]) % 2
            # return probability of op
            return op_to_pr[pt.bsf_to_pauli(op)]

        @functools.lru_cache(256)  # 256=2^8 is enough for fixed prob_dist
        def create_q_node(self, prob_dist, fs, compass_direction=None):
            """
            Return qubit tensor.

            :param prob_dist: Tuple of probability distribution in the format (P(I), P(X), P(Y), P(Z)).
            :type prob_dist: 4-tuple of float
            :param fs: Two qubit Paulis, e.g. ('I', None) or ('Z', 'I').
            :type fs: 2-tuple of str
            :param compass_direction: Location of tensor relative to squared lattice (figure 3), e.g. 'n', 'se' or None.
            :type compass_direction: str
            :return: Qubit tensor
            :rtype: numpy.array (4d)
            """

            def _node_shapes(direction=None):
                """Return shape of upper and lower tensors including dummy indices."""
                return {  # direction order n,e,s,w
                    'n': ((1, 4, 1, 4), (1, 4, 4, 4)),
                    'ne': (None, (1, 1, 4, 4)),
                    'e': ((1, 1, 1, 4), None),
                    'se': ((4, 1, 1, 4), None),
                    's': ((4, 4, 1, 4), (1, 4, 1, 4)),
                    'sw': ((4, 1, 1, 1), None),
                    'w': ((4, 4, 1, 1), (1, 4, 4, 1)),
                    'nw': ((1, 4, 1, 1), (1, 4, 4, 1)),
                }.get(direction, ((4, 4, 1, 4), (1, 4, 4, 4)))

            # create tensors with values
            nodes = []
            for shape, f, permute in zip(_node_shapes(compass_direction), fs, (False, True)):
                # stabilizers on lower qubits are permutable
                assert (shape is None) == (f is None), 'Restricted Paulis do not match shapes.'
                if shape is None:  # add dummy tensor
                    nodes.append(np.ones((1, 1, 1, 1), dtype=np.float64))
                else:  # add qubit tensor
                    node = np.empty(shape, dtype=np.float64)
                    for n, e, s, w in np.ndindex(node.shape):
                        node[(n, e, s, w)] = self.q_node_value(prob_dist, f, n, e, s, w, permute)
                    nodes.append(node)
            # merge upper and lower tensors
            node = np.einsum('nesw,sESW->neESwW', nodes[0], nodes[1]).reshape(
                (
                    nodes[0].shape[0],  # n
                    nodes[0].shape[1] * nodes[1].shape[1],  # eE
                    nodes[1].shape[2],  # S
                    nodes[0].shape[3] * nodes[1].shape[3]  # wW
                )
            )
            # multiply dimension-16 legs with deltas to reduce dimension to 4
            if node.shape[1] == 16:
                node = np.einsum('nesw,Ee->nEsw', node, tt.tsr.delta((4, 4, 4)).reshape((4, 16)))
            if node.shape[3] == 16:
                node = np.einsum('nesw,Ww->nesW', node, tt.tsr.delta((4, 4, 4)).reshape((4, 16)))
            return node
