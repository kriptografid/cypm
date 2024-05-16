# CyPM - Python module to analyze Column Parity Mixers (CPM)

Column Parity Mixers (CPM) are the linear functions used to provide diffusion in the symmetric-key primitives,
specifically permutation-based cryptography such as Keccak (SHA-3), Xoodoo, and Ascon. The description of CPM
is available in the following paper:

- Ko Stoffelen and Joan Daemen. [Column Parity Mixers](https://doi.org/10.13154/tosc.v2018.i1.126-159). IACR Transactions on Symmetric Cryptology (ToSC).
Volume 2018. Issue 1.


## How to Use

This module runs on top of SageMath. To use it, simply clone or download this repository and
run the following from the sage terminal

    sage: load("cypm.py")

Note that for the above command, you must run the sage terminal from the same folder as this repository. Otherwise, change `"cypm.py"`
to its absolute path location.


## Examples

The CPM is defined by a positive integer `m` and a squared binary matrix `Z` called *parity folding matrix*.

    sage: load("cypm.py")
    sage: Z = [
        [1,0,0,0,0,0,0,1],
        [1,1,0,0,0,0,0,0],
        [0,1,1,0,0,0,0,0],
        [0,0,1,1,0,0,0,0],
        [0,0,0,1,1,0,0,0],
        [0,0,0,0,1,1,0,0],
        [0,0,0,0,0,1,1,0],
        [0,0,0,0,0,0,1,1],
    ]
    sage: theta = CPM(5, Z)
    sage: theta
    Column Parity Mixer with 5 rows and 8 columns

The properties of `theta` can be inquired easily, such as checking is invertibility and compute its inverse.

    sage: theta.is_invertible()
    True
    sage: thetainv = theta.inverse()  # compute the inverse of theta
    sage: thetainv
    Column Parity Mixer with 5 rows and 8 columns
    sage: thetainv.parity_folding_matrix()  # the parity folding matrix of the inverse of theta
    [1 1 0 0 0 0 0 0]
    [0 1 1 0 0 0 0 0]
    [0 0 1 1 0 0 0 0]
    [0 0 0 1 1 0 0 0]
    [0 0 0 0 1 1 0 0]
    [0 0 0 0 0 1 1 0]
    [0 0 0 0 0 0 1 1]
    [1 0 0 0 0 0 0 1]

One can check the `theta` effect on certain input matrices,

    sage: A = [
        [0,0,0,0,0,1,0,0],
        [1,0,0,0,0,1,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,1,0,0,1,0,0],
        [0,1,1,0,0,1,0,0],
    ] 
    sage: theta.effect(A)
    (0, 1, 0, 0, 1, 1, 0, 1)
    sage: theta.affected_columns(A)
    [1, 4, 5, 7]
    sage: theta(A)  # compute the image of A under theta
    [0 1 0 0 1 0 0 1]
    [1 1 0 0 1 0 0 1]
    [0 1 0 0 1 0 0 1]
    [0 1 1 0 1 0 0 1]
    [0 0 1 0 1 0 0 1]

or compute the composition of `theta` with other CPM,

    sage: Z = [
        [0, 1, 1, 0, 0, 0, 1, 1],
        [1, 1, 0, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 1],
        [0, 1, 0, 1, 1, 0, 0, 1],
        [0, 1, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 1, 0, 1, 1, 0],
        [1, 0, 1, 1, 0, 0, 0, 1],
        [0, 1, 0, 1, 0, 1, 1, 1]
    ]
    sage: gamma = CPM(5, Z)
    sage: gamma
    sage: alpha = theta.compose(gamma)
    sage: alpha.parity_folding_matrix()
    [1 1 0 1 0 1 1 0]
    [1 0 1 0 0 0 1 1]
    [1 0 1 0 0 1 1 1]
    [0 0 1 1 0 1 0 1]
    [0 1 0 0 0 0 0 1]
    [0 1 1 1 0 0 0 0]
    [0 0 0 1 0 0 0 0]
    [1 0 1 1 0 0 1 0]
    sage: theta_3 = theta^3  # composition of theta with itself 3-times
    sage: theta.compose(theta).compose(theta) == theta_3
    True

as well as computing the order of the CPMs

    sage: theta.order()
    8
    sage: theta_3.order()
    8
    sage: gamma.order()
    85


## License

CyPM is released under GNU General Public License version 3 (GPLv3)


## Copyright

[PT Kriptograf Indo Teknologi](https://kriptograf.id) - [contact@kriptograf.id](mailto:contact@kriptograf.id)
