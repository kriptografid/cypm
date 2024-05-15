"""
CyPM - SageMath module to analyze Column Parity Mixers (CPM)

Copyright (C) 2024 -  PT Kriptograf Indo Teknologi <contact@kriptograf.id>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from sage.all import GF, vector, Matrix
from sage.matrix.special import identity_matrix
from sage.misc.functional import is_even, is_odd


class CPM(object):
    def __init__(self, m, Z):
        """
        Construct an instance of Column Parity Mixer (CPM)

        INPUT:

        - ``m`` -- a positive integer
        - ``Z`` -- parity folding matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [            \
                [1,0,0,0,0,0,0,1], \
                [1,1,0,0,0,0,0,0], \
                [0,1,1,0,0,0,0,0], \
                [0,0,1,1,0,0,0,0], \
                [0,0,0,1,1,0,0,0], \
                [0,0,0,0,1,1,0,0], \
                [0,0,0,0,0,1,1,0], \
                [0,0,0,0,0,0,1,1], \
            ]
            sage: theta = CPM(5, Z)
            sage: theta
            Column Parity Mixer with 5 rows and 8 columns
            sage: theta.is_invertible()
            True
            sage: theta.order()
            8
            sage: theta.associated_matrix()
            [0 0 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]
            [0 0 0 0 1 0 0 0]
            [0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 1 0]
        """
        if m < 1:
            raise ValueError("m must be >= 1")

        if isinstance(Z, (list, tuple)):
            Z = Matrix(GF(2), Z)

        if Z.base_ring() != GF(2):
            raise ValueError("parity folding matrix must be defined over GF(2)")

        if not Z.is_square():
            raise TypeError("Z must be a square matrix")

        self._nrows = m
        self._ncols = Z.ncols()
        self._Z = Z

    @staticmethod
    def column_parity(A):
        """
        Return the column parity of the matrix A

        INPUT:

        - ``A`` -- a Matrix over GF(2)

        EXAMPLES::

            sage: load("cypm.py")
            sage: A = [             \
                [0,0,0,0,0,1,0,0],  \
                [1,0,0,0,0,1,0,0],  \
                [0,0,0,0,0,1,0,0],  \
                [0,0,1,0,0,1,0,0],  \
                [0,1,1,0,0,1,0,0],  \
            ]
            sage: CPM.column_parity(A)
            (1, 1, 0, 0, 0, 1, 0, 0)
        """
        A = Matrix(GF(2), A)
        return sum(A)  # sum all rows of A

    @staticmethod
    def expanded_column_parity(A):
        """
        Return the expanded column parity of the matrix A

        INPUT:

        - ``A`` -- a Matrix over GF(2)

        EXAMPLES::

            sage: load("cypm.py")
            sage: A = [             \
                [0,0,0,0,0,1,0,0],  \
                [1,0,0,0,0,1,0,0],  \
                [0,0,0,0,0,1,0,0],  \
                [0,0,1,0,0,1,0,0],  \
                [0,1,1,0,0,1,0,0],  \
            ]
            sage: CPM.expanded_column_parity(A)
            [1 1 0 0 0 1 0 0]
            [1 1 0 0 0 1 0 0]
            [1 1 0 0 0 1 0 0]
            [1 1 0 0 0 1 0 0]
            [1 1 0 0 0 1 0 0]
        """
        A = Matrix(GF(2), A)
        m = A.nrows()
        row = CPM.column_parity(A)
        return Matrix(GF(2), [row for _ in range(m)])

    @staticmethod
    def odd_columns(A):
        """
        Return the odd columns of A

        INPUT:

        - ``A`` -- a Matrix over GF(2)

        EXAMPLES::

            sage: load("cypm.py")
            sage: A = [             \
                [0,0,0,0,0,1,0,0],  \
                [1,0,0,0,0,1,0,0],  \
                [0,0,0,0,0,1,0,0],  \
                [0,0,1,0,0,1,0,0],  \
                [0,1,1,0,0,1,0,0],  \
            ]
            sage: CPM.odd_columns(A)
            [0, 1, 5]
        """
        col_parity = CPM.column_parity(A)
        return col_parity.nonzero_positions()

    @staticmethod
    def even_columns(A):
        """
        Return the even columns of A

        INPUT:

        - ``A`` -- a Matrix over GF(2)

        EXAMPLES::

            sage: load("cypm.py")
            sage: A = [             \
                [0,0,0,0,0,1,0,0],  \
                [1,0,0,0,0,1,0,0],  \
                [0,0,0,0,0,1,0,0],  \
                [0,0,1,0,0,1,0,0],  \
                [0,1,1,0,0,1,0,0],  \
            ]
            sage: CPM.even_columns(A)
            [2, 3, 4, 6, 7]
        """
        A = Matrix(GF(2), A)
        return [e for e in range(A.ncols()) if e not in CPM.odd_columns(A)]

    def nrows(self):
        """
        Return the no. of rows of the CPM

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: theta.nrows()
            4
        """
        return self._nrows

    def ncols(self):
        """
        Return the no. of columsn of the CPM

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: theta.ncols()
            3
        """
        return self._ncols

    def parity_folding_matrix(self):
        """
        Return the parity folding matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: theta.parity_folding_matrix()
            [1 1 0]
            [0 1 1]
            [1 0 1]
        """
        return self._Z

    def effect(self, A):
        """
        Return the effect of `A`

        INPUT:

        - ``A`` -- input matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: A = [[1, 1, 1], [0, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: theta.effect(A)
            (1, 1, 0)
        """
        Z = self.parity_folding_matrix()
        return CPM.column_parity(A) * Z

    def expanded_effect(self, A):
        """
        Return the expanded effect of `A`

        INPUT:

        - ``A`` -- input matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: A = [[1, 1, 1], [0, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: theta.expanded_effect(A)
            [1 1 0]
            [1 1 0]
            [1 1 0]
            [1 1 0]
        """
        Z = self.parity_folding_matrix()
        return CPM.expanded_column_parity(A) * Z

    def affected_columns(self, A):
        """
        Return a list of affected columns

        INPUT:

        - ``A`` -- input matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: A = [[1, 1, 1], [0, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: theta.affected_columns(A)
            [0, 1]
        """
        return self.effect(A).nonzero_positions()

    def unaffected_columns(self, A):
        """
        Return a list of unaffected columns

        INPUT:

        - ``A`` -- input matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
            sage: theta = CPM(4, Z)
            sage: A = [[1, 1, 1], [0, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: theta.unaffected_columns(A)
            [2]
        """
        return [i for i in range(self.ncols()) if i not in self.affected_columns(A)]

    def __call__(self, A):
        """

        INPUT:

        - ``A`` -- input matrix

        TESTS::

            sage: load("cypm.py")
            sage: Z = [            \
                [1,0,0,0,0,0,0,1], \
                [1,1,0,0,0,0,0,0], \
                [0,1,1,0,0,0,0,0], \
                [0,0,1,1,0,0,0,0], \
                [0,0,0,1,1,0,0,0], \
                [0,0,0,0,1,1,0,0], \
                [0,0,0,0,0,1,1,0], \
                [0,0,0,0,0,0,1,1], \
            ]
            sage: A = [            \
                [0,0,0,0,0,1,0,0], \
                [1,0,0,0,0,1,0,0], \
                [0,0,0,0,0,1,0,0], \
                [0,0,1,0,0,1,0,0], \
                [0,1,1,0,0,1,0,0], \
            ]
            sage: theta = CPM(5, Z)
            sage: theta
            Column Parity Mixer with 5 rows and 8 columns
            sage: theta(A)
            [0 1 0 0 1 0 0 1]
            [1 1 0 0 1 0 0 1]
            [0 1 0 0 1 0 0 1]
            [0 1 1 0 1 0 0 1]
            [0 0 1 0 1 0 0 1]
        """
        A = Matrix(GF(2), A)
        return A + self.expanded_effect(A)

    def compose(self, other):
        """
        Return the composition of this CPM with the `other` CPM

        INPUT:

        - ``other`` -- an instance of CPM

        EXAMPLES::

            sage: load("cypm.py")
            sage: n = 7
            sage: theta_p = CPM(4, random_matrix(GF(2), n, n))
            sage: theta = CPM(4, random_matrix(GF(2), n, n))
            sage: alpha = theta_p.compose(theta)
            sage: A = random_matrix(GF(2), alpha.nrows(), alpha.ncols())
            sage: alpha(A) == theta(theta_p(A))
            True

        TESTS::

            sage: theta_p = CPM(5, random_matrix(GF(2), n, n))
            sage: theta = CPM(5, random_matrix(GF(2), n, n))
            sage: alpha = theta_p.compose(theta)
            sage: A = random_matrix(GF(2), alpha.nrows(), alpha.ncols())
            sage: alpha(A) == theta(theta_p(A))
            True
        """
        if not isinstance(other, CPM):
            raise TypeError("other must be an instance of CPM")

        if self.nrows() != other.nrows() or self.ncols() != other.ncols():
            raise ValueError("the dimension of the CPM must be equal")

        m = self.nrows()
        if is_even(m):
            Z = self.parity_folding_matrix() + other.parity_folding_matrix()
        else:
            I = identity_matrix(GF(2), self.ncols())
            Z = ((self.parity_folding_matrix() + I) * (other.parity_folding_matrix() + I)) + I

        return CPM(m, Z)

    def order(self):
        """
        Return the order of the CPM

        EXAMPLES::

            sage: load("cypm.py")
            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: theta.order()
            7
        """
        m = self.nrows()
        if is_even(m):
            return 2
        return self.associated_matrix().multiplicative_order()

    def associated_matrix(self):
        """
        Return the associated matrix, i.e. the sum of parity folding matrix and the identity matrix

        EXAMPLES::

            sage: load("cypm.py")
            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: theta.associated_matrix()
            [0 1 0]
            [1 1 1]
            [0 1 1]
        """
        return self.parity_folding_matrix() + identity_matrix(GF(2), self.ncols())

    def is_invertible(self):
        """
        Return `True` if this CPM is invertible

        EXAMPLES::

            sage: load("cypm.py")
            sage: c = CPM(2, [[1, 0, 1], [1, 1, 1], [0, 1, 0]])
            sage: c.is_invertible()
            True
            sage: Z = [            \
                [1,0,0,0,0,0,0,1], \
                [1,1,0,0,0,0,0,0], \
                [0,1,1,0,0,0,0,0], \
                [0,0,1,1,0,0,0,0], \
                [0,0,0,1,1,0,0,0], \
                [0,0,0,0,1,1,0,0], \
                [0,0,0,0,0,1,1,0], \
                [0,0,0,0,0,0,1,1], \
            ]
            sage: theta = CPM(5, Z)
            sage: theta.is_invertible()
            True
        """
        m = self.nrows()
        return is_even(m) or not self.associated_matrix().is_singular()

    def inverse(self):
        """
        Return the inverse of this CPM

        EXAMPLES::

            sage: load("cypm.py")
            sage: Z = [            \
                [1,0,0,0,0,0,0,1], \
                [1,1,0,0,0,0,0,0], \
                [0,1,1,0,0,0,0,0], \
                [0,0,1,1,0,0,0,0], \
                [0,0,0,1,1,0,0,0], \
                [0,0,0,0,1,1,0,0], \
                [0,0,0,0,0,1,1,0], \
                [0,0,0,0,0,0,1,1], \
            ]
            sage: theta = CPM(5, Z)
            sage: thetainv = theta.inverse()
            sage: I = theta.compose(thetainv)
            sage: I.associated_matrix()
            [1 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]
            [0 0 0 0 1 0 0 0]
            [0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 1]
        """
        if not self.is_invertible():
            raise ValueError("CPM is not invertible")

        m = self.nrows()
        if is_even(m):
            return self

        Z = self.associated_matrix().inverse() - identity_matrix(GF(2), self.ncols())
        return CPM(m, Z)

    def is_equal(self, other):
        """
        Return `True` if this CPM is equal to the `other` CPM

        INPUT:

        - ``other`` -- an instance of CPM

        EXAMPLES::

            sage: load("cypm.py")
            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: thetap = CPM(5, [[1, 1, 1], [1, 0, 1], [0, 1, 0]])
            sage: theta.is_equal(thetap)
            False
        """
        if not isinstance(other, CPM):
            raise TypeError("the other input must be an instance of CPM")
        return self.nrows() == other.nrows() and self.parity_folding_matrix() == other.parity_folding_matrix()

    def is_identity(self):
        """
        Return `True` if this CPM is an identity map

        EXAMPLES::

            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: theta.is_identity()
            False
            sage: theta.order()
            7
            sage: (theta**7).is_identity()
            True
        """
        return self.parity_folding_matrix().is_zero()

    def is_involution(self):
        """
        Return `True` if this CPM is an involution

        EXAMPLES::

            sage: load("cypm.py")
            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: theta.is_involution()
            False
            sage: thetap = CPM(4, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: thetap.is_involution()
            True
            sage: thetap.compose(thetap).is_identity()
            True
        """
        m = self.nrows()
        if is_even(m):
            return True
        return self.compose(self).is_identity()

    def __repr__(self):
        return f"Column Parity Mixer with {self.nrows()} rows and {self.ncols()} columns"

    def __eq__(self, other):
        """
        Overload equality comparison operator

        EXAMPLES::

            sage: load("cypm.py")
            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: thetap = CPM(5, [[1, 1, 1], [1, 0, 1], [0, 1, 0]])
            sage: theta == thetap
            False
            sage: theta == theta
            True
        """
        return self.is_equal(other)

    def __pow__(self, e):
        """
        Return the iterative composition of this CPM for `e` times

        INPUT:

        - ``e`` -- a positive integer

        EXAMPLES::

            sage: theta = CPM(5, [[1, 1, 0], [1, 0, 1], [0, 1, 0]])
            sage: theta.order()
            7
            sage: theta**8 == theta
            True
        """
        if e < 1:
            raise ValueError("e must be >= 1")

        m, Z = self.nrows(), self.parity_folding_matrix()
        theta = CPM(m, Z)
        for _ in range(e-1):
            theta = theta.compose(self)
        return theta
