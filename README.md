# C++ Lab Assignment: Linear Algebra with Eigen Library

## Objective
In this discussion, we will review some 
fundamental linear algebra concepts and how to use the Eigen library in C++.

Linear algebra is important in many areas of science, and in particular, quantum chemistry
is based almost entirely on some linear algebra operations (namely, tensor contractions and
matrix diagonalization).

## Background and Review

### Vectors
- **Dot Product**: Computes the projection of one vector onto another, resulting in a scalar value.
  Is the sum of products of two vectors (which must be the same length).
  $$
  \vec{u} \cdot \vec{v} = \sum_{i=1}^{n} u_i v_i
  $$
- **Outer Product**: Produces a matrix by multiplying a column vector with a row vector. Can be thought
  of as basically considering the two vectors to be $N \times 1$ and a $1 \times N$ matrices, and doing
  matrix multiplication.
  
  $$
    (\vec{u} \otimes \vec{v})_{ij} = u_iv_j
  $$

### Matrices

Matrices have two indices, here row $i$ and column $j$

- **Addition**: Add two matrices element-wise.

  $$
    (\mathbf{a} + \mathbf{b})_{ij} = a_{ij}b_{ij}
  $$

- **Multiplication**: Multiplies rows of the left matrix with columns of the right. Can think of each
  each element of the resulting matrix to be a vector dot product between a row of $A$ and a column of $B$.

  $$
  (AB)_{ij} = \sum_k A_{ik}B_{kj}
  $$

- **Transpose**: Flip a matrix across its diagonal

  $$
  (A^T)_{ij} = A_{ji}
  $$

- **Conjugate transpose**: Similar to transpose, but for complex matrices.

  The *conjugate* of a complex number has the same real part. The imaginary part is
  the same magnitude, but opposite sign. For example, the conjugate of $3.0+4.1i$
  is $3.0-4.1i$.

  The conjugate transpose of a matrix has symmetric real parts, but conjugate imaginary parts.

  In quantum chemistry and physics, conjugate transpose is usually denoted with a dagger $\dagger$
  rather than a $T$ for transpose.

  $$
  (A^\dagger)_{ij} = \bar{A_{ji}}
  $$

  where the overbar represents conjugation.

- **Diagonal matrix**: A matrix $\mathbf{A}$ is diagonal if off-diagonal elements are zero.
  $$
  A_{i \neq j} = 0
  $$

  or, more visually, a diagonal matrix

  $$
    D =
    \begin{bmatrix}
    d_{11} & 0 & 0 & 0 \\
    0 & d_{22} & 0 & 0 \\
    0 & 0 & d_{33} & 0 \\
    0 & 0 & 0 & d_{44}
    \end{bmatrix}
    $$

    and a non-diagonal matrix
    $$
    A =
    \begin{bmatrix}
    a_{11} & a_{12} & a_{13} & a_{14} \\
    a_{21} & a_{22} & a_{23} & a_{24} \\
    a_{31} & a_{32} & a_{33} & a_{34} \\
    a_{41} & a_{42} & a_{43} & a_{44}
    \end{bmatrix}
    $$

### Symmetric and Hermitian Matrices

A matrix is symmetric if it is equal to its own transpose. That is,
$\mathbf{A}_{ij} = \mathbf{A}_{ji}$ Below is an example of a symmetric matrix.

  $$
    S =
    \begin{bmatrix}
    2.0 & 1.4 & 8.2 & 0.0 \\
    1.4 & 3.0 & 0.9 & 7.5 \\
    8.2 & 0.9 & 4.0 & 9.8 \\
    0.0 & 7.5 & 9.8 & 5.0
    \end{bmatrix}
  $$

Typically, when referring to a "symmetric" matrix, we are referring to a *real* (not complex)
symmetric matrix.

The extension of this for complex matrices is a *self-adjoint* or *hermitian* matrix.
A hermitian matrix is equal to its own *conjugate transpose*.

$$
A =
\begin{bmatrix}
1.2 + 0.0i & 2.3 + 1.1i & 3.4 + 0.5i & 4.1 + 0.3i \\
2.3 + 1.1i & 5.6 + 0.0i & 6.7 + 1.2i & 7.8 + 0.4i \\
3.4 + 0.5i & 6.7 + 1.2i & 8.9 + 0.0i & 9.1 + 1.3i \\
4.1 + 0.3i & 7.8 + 0.4i & 9.1 + 1.3i & 2.2 + 0.0i
\end{bmatrix}
$$

Symmetric and hermitian matrices have special properties, which we will see below.


### Matrix Diagonalization

A very important operation that is done on matrices is diagonalizing them. That is, we decompose a single
matrix into a diagonal matrix and other matrices.

$$
  A = \mathbf{P} \mathbf{D} \mathbf{P}^{-1}
$$

where $\mathbf{P}^{-1}$ is the inverse of matrix $\mathbf{P}$ (that is, $\mathbf{PP}^{-1} = \mathbf{I}$, the
identity matrix, which is a diagonal matrix with all ones on the diagonal).

Importantly, the values along the diagonal elements of diagonal matrix $\mathbf{D}$ are called *eigenvalues*
and the columns of $\mathbf{P}$ are called *eigenvectors*. Eigenvalues and eigenvectors of various matrices
have enormous importance in quantum chemistry, as well as many areas of science.

One particular thing to note is that eigenvalues and eigenvectors of matrices are often complex. However,
**eigenvalues and eigenvectors of hermitian matrices are real**.

## Some conceptual questions

1. Is a real, symmetric matrix also hermitian?
2. What seems special about the diagonal of the hermitian matrix given above? Is this always true? Why or why not?
3. Why, as a programmer, would we care whether the eigenvalues/eigenvectors of a matrix are real or complex?

## Programming Assignment

For this discussion, we will be using the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) package.
One convenient feature is that Eigen is header-only, so does not require compiling a library.

Create a directory for this discussion. Inside that directory, download this file and unpack it: https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2

Open the `main.cpp` source file, and then lets have some fun with linear algebra.

### Some Eigen notes

Some common matrix and vector types (you can see a pattern):

* `MatrixXd` - arbitrary-sized double precision
* `Matrix3d` - $3\times3$ double precision
* `Vector3d` - double precision vector of length 3
* `MatrixXcd` - arbitrary-sized complex double precision

There's lots of functions associated with Eigen. Of particular use in this discussion
is the `.transpose()` method of matrices.

If you have coordinates already, it's easy to construct an Eigen matrix using brace initialization:

```c++
    Eigen::MatrixXd m{
       {  1.0,  2.0,  3.0,  4.0},
       {  5.0,  6.0,  7.0,  8.0},
       {  9.0, 10.0, 11.0, 12.0},
       { 13.0, 14.0, 15.0, 16.0},
    };
```


### Rotating and Modifying Coordinates

A 3D rotation can be represented by a $3 \times 3$ rotation matrix $\mathbf{R}$.
An $N \times 3$ matrix of points $\mathbf{M}$ can be rotated by multplication with this
rotation matrix.

  $$
  \mathbf{M}' = \mathbf{M} \mathbf{R}^T
  $$


For example, a rotation matrix that rotates the coordinates by 90 degrees clockwise around the z-axis would be

  $$
  \begin{bmatrix}
   0.0  &  1.0  & 0.0\\
  -1.0  &  0.0 &  0.0\\
   0.0  &  0.0 &  1.0
  \end{bmatrix}
  $$

Convince yourself graphically that this is true.

In the `main.cpp`, I have given you the coordinates of two molecules - water, and formaldehyde.

Rotate the water molecule by 90 cdegrees using a rotation matrix.

Create a rotation matrix for rotating by 90 degrees counterclockwise, and apply that to the already-rotated
coordinates. Verify you end up back where you started.

### Diagonalizing a matrix

The Eigen library includes classes for diagonalizing matrices. In particular, there is the
[`SelfAdjointEigenSolver`](https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html#af9cf17478ced5a7d5b8391bb10873fac) for solving self-adjoint/hermitian matrices,
and the more general [`EigenSolver`](https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html#a0ccaeb4f7d44c18af60a7b3a1dd91f7a).
Both of these classes are templated (that is, you give it the matrix type). I have linked above to the documentation
that shows examples of using these classes. In general, you construct the class with your matrix to be diagonalized,
then call `.eigenvalues()` and `.eigenvectors()`.

**Question** - If the matrix type for the general `EigenSolver` is `MatrixXd`, what do you expect the return type
of the eigenvalues and eigenvectors to be?

**Question** - How does the example code for the `SelfAdjointEigenSolver` create a symmetric matrix?

In `main.cpp`, I have given you a square matrix. Determine which eigensolver to use and form the
eigenvectors and eigenvalues. Show that the decomposition is true by forming $\mathbf{P} \mathbf{D} \mathbf{P}^{-1}$.
You should end up with your original matrix again.

### Calculating the frequencies of a molecule

One area (among many) where diagonalization shows up is in the calculation of vibrational frequencies.
These frequencies are part of a molecule's infrared (IR) spectrum. Here, for example, is an IR
spectrum of formaldehyde: https://webbook.nist.gov/cgi/cbook.cgi?ID=C50000&Type=IR-SPEC&Index=0#IR-SPEC

To compute this with quantum chemistry, we need software that will compute the *hessian* matrix.
This is a square matrix that represents
the second derivative of energy of the molecule with respect to the atomic positions.

  $$
  \begin{bmatrix}
   F_{x1,x1}   &  F_{x1,y1}  &  F_{x1,z1}  &  F_{x1,x2}  &  \cdots  &  F_{x1,zN}\\
   F_{y1,x1}   &  F_{y1,y1}  &  F_{y1,z1}  &  F_{y1,x2}  &  \cdots  &  F_{y1,zN}\\
   F_{z1,x1}   &  F_{z1,y1}  &  F_{z1,z1}  &  F_{z1,x2}  &  \cdots  &  F_{z1,zN}\\
   F_{x2,x1}   &  F_{x2,y1}  &  F_{x2,z1}  &  F_{x2,x2}  &  \cdots  &  F_{x2,zN}\\
   \vdots \\
   F_{zN,x1}   &  F_{zN,y1}  &  F_{zN,z1}  &  F_{zN,x2}  &  \cdots  &  F_{zN,zN}\\
  \end{bmatrix}
  $$

where $N$ is the number of atoms.

We will need to diagonalize this matrix, but first, we need to write some other functions.
Normally, this matrix is computed by quantum chemistry software, but I have provided you
some hessians in some files. So we need to write a function to read them in.

#### Reading the hessian from a file

Write a function that reads in the hessian matrix from a file. The funtion needs to take two arguments - the path
to a file, and the number of atoms. It should return an Eigen matrix of some kind. What is the size of that matrix
in terms of the number of atoms?

There are two files in the discussion repo - `water_hessian.dat` and `formaldehyde_hessian.dat`. Read both of these
into matrices and verify (by printing them out) that they were read in properly. **Hint** - Eigen matrices directly
support printing (ie, `std::cout << mat << std::end;`)

#### Mass-weighting the hessian

Next, we must modify the hessian.

First, you should create a vector (`std::vector`) of the atomic masses for each molecule.
I have provided a lookup map already. The ordering of the atoms in each molecule is given in a comment above the coordinates.

The next part is a little tricky. You must scale each element of the hessian matrix by the sqrt of the product of masses.
This is the *mass-weighted* hessian.

$$
H_{ij}^{mw} = \frac{H_{ij}}{\sqrt{m_im_j}}
$$

**Note that the hessian matrix contains x, y, and z coordinates**. That is, $x_1$, $y_1$, and $z_1$ should all be scaled
by the mass of the first atom.

**Question** - Given index $i$ in the hessian matrix, what is a simple expression for the atom number?
If $i=0$, then it's atom $0$. If $i=1$, it's atom $0$. If $i=3$, it's atom $1$, etc.

With that, you can loop over the elements of the hessian and form a new mass-weighted hessian matrix.


#### Diagonalizing the mass-weighted hessian

Next, you need to diagonalize this matrix. Is this matrix symmetric?

For now, we only care about the eigenvalues. The eigenvalues represent the frequency of the vibration. This is
directly related to it's infrared spectrum.

Once you have the eigenvalues, you need to convert the values to wavenumbers (in $\mathrm{cm}^{-1}$). To do this,
use the following formula:

$$
\tilde{v}_i = 5140.484532 * \sqrt{\lambda_i}
$$

where $\lambda_i$ is the corresponding eigenvalue.

#### Comparison

Compute the vibrational frequencies for formaldehyde and compare them with the spectrum [here](https://webbook.nist.gov/cgi/cbook.cgi?ID=C50000&Type=IR-SPEC&Index=0#IR-SPEC). How well does it compare?
