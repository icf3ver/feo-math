use crate::{F32Fmt, SignOps, /*Sum,*/ Two};

use {
    super::{
        SqMatrix,
        vector4::Vector4,
    },
    crate::{
        Construct,
        One, Zero
    },
    std::{
        fmt,
        ops::{Add, Mul, Sub, Div, Rem, Neg}
    },
};

/// A 4x4 matrix with row major order.
#[derive(PartialEq, Clone, Copy)]
pub struct Matrix4<T>{pub m:[[T; 4]; 4]}

impl<T> Matrix4<T> {
    /// A cleaner way of making a new Matrix4.
    /// Returns a new Matrix4<T>
    /// # Arguments
    /// * `row0` - A 4 element array representing the first row of the matrix
    /// * `row2` - A 4 element array representing the second row of the matrix
    /// * `row3` - A 4 element array representing the third row of the matrix
    /// * `row4` - A 4 element array representing the fourth row of the matrix
    pub fn new(
            row0: [T; 4], 
            row1: [T; 4], 
            row2: [T; 4], 
            row3: [T; 4]) -> Self {
        Matrix4 {
            m: [
                row0,
                row1,
                row2,
                row3
            ]
        }
    }
    
    #[allow(clippy::too_many_arguments)]
    pub fn new_val(
        r0c0: T, r0c1: T, r0c2: T, r0c3: T,
        r1c0: T, r1c1: T, r1c2: T, r1c3: T,
        r2c0: T, r2c1: T, r2c2: T, r2c3: T,
        r3c0: T, r3c1: T, r3c2: T, r3c3: T) -> Self {
        Matrix4 {
            m: [
                [r0c0, r0c1, r0c2, r0c3],
                [r1c0, r1c1, r1c2, r1c3],
                [r2c0, r2c1, r2c2, r2c3],
                [r3c0, r3c1, r3c2, r3c3]
            ]
        }
    }

    /// one: T equivalent of one zero: T equivalent of zero
    pub fn identity() -> Self
    where T: Zero + One {
        Self::ONE
    }

    #[allow(dead_code)]
    pub fn determinant(&self) -> T 
    where T: Add<T, Output = T> + Mul<T, Output = T> + Sub<T, Output = T> + Copy {
        self.m[0][0] * (
                (self.m[1][1] * self.m[2][2] * self.m[3][3]) +
                (self.m[1][2] * self.m[2][3] * self.m[3][1]) +
                (self.m[1][3] * self.m[2][1] * self.m[3][2])
                - (self.m[1][3] * self.m[2][2] * self.m[3][1])
                - (self.m[1][2] * self.m[2][1] * self.m[3][3])
                - (self.m[1][1] * self.m[2][3] * self.m[3][2])
            ) - self.m[1][0] * (
                (self.m[0][1] * self.m[2][2] * self.m[3][3]) +
                (self.m[0][2] * self.m[2][3] * self.m[3][1]) +
                (self.m[0][3] * self.m[2][1] * self.m[3][2])
                - (self.m[0][3] * self.m[2][2] * self.m[3][1])
                - (self.m[0][2] * self.m[2][1] * self.m[3][3])
                - (self.m[0][1] * self.m[2][3] * self.m[3][2])
            ) + self.m[2][0] * (
                (self.m[0][1] * self.m[1][2] * self.m[3][3]) +
                (self.m[0][2] * self.m[1][3] * self.m[3][1]) +
                (self.m[0][3] * self.m[1][1] * self.m[3][2])
                - (self.m[0][3] * self.m[1][2] * self.m[3][1])
                - (self.m[0][2] * self.m[1][1] * self.m[3][3])
                - (self.m[0][1] * self.m[1][3] * self.m[3][2])
            ) - self.m[3][0] * (
                (self.m[0][1] * self.m[1][2] * self.m[2][3]) +
                (self.m[0][2] * self.m[1][3] * self.m[2][1]) +
                (self.m[0][3] * self.m[1][1] * self.m[2][2])
                - (self.m[0][3] * self.m[1][2] * self.m[2][1])
                - (self.m[0][2] * self.m[1][1] * self.m[2][3])
                - (self.m[0][1] * self.m[1][3] * self.m[2][2])
            )
    }

    // transpose
    // matrix (a, b, c, d)
    //        (e, f, g, h)
    //        (i, j, k, l)
    //        (m, n, o, p)
    //
    // matrix (a, d, g)
    //        (b, e, h)
    //        (c, f, i) // nd
    pub fn transpose(&self) -> Matrix4<T>
    where T: Copy {
        Matrix4 {
            m: [
                [self.m[0][0], self.m[1][0], self.m[2][0], self.m[3][0]],
                [self.m[0][1], self.m[1][1], self.m[2][1], self.m[3][1]],
                [self.m[0][2], self.m[1][2], self.m[2][2], self.m[3][2]],
                [self.m[0][3], self.m[1][3], self.m[2][3], self.m[3][3]]
            ]
        }
    }

    #[allow(clippy::possible_missing_comma)]
    pub fn cofactor(self) -> Self
    where T: Mul<T, Output = T> + Neg<Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Copy {
        Matrix4 {
            m: [
                [
                    (self.m[1][1] * self.m[2][2] * self.m[3][3]) +
                    (self.m[2][1] * self.m[3][2] * self.m[1][3]) +
                    (self.m[3][1] * self.m[1][2] * self.m[2][3])
                    - (self.m[3][1] * self.m[2][2] * self.m[1][3])
                    - (self.m[2][1] * self.m[1][2] * self.m[3][3])
                    - (self.m[1][1] * self.m[3][2] * self.m[2][3]),

                    - (self.m[1][0] * self.m[2][2] * self.m[3][3])
                    - (self.m[2][0] * self.m[3][2] * self.m[1][3])
                    - (self.m[3][0] * self.m[1][2] * self.m[2][3]) + 
                    (self.m[3][0] * self.m[2][2] * self.m[1][3]) +
                    (self.m[2][0] * self.m[1][2] * self.m[3][3]) +
                    (self.m[1][0] * self.m[3][2] * self.m[2][3]),

                    (self.m[1][0] * self.m[2][1] * self.m[3][3]) + 
                    (self.m[2][0] * self.m[3][1] * self.m[1][3]) + 
                    (self.m[3][0] * self.m[1][1] * self.m[2][3])
                    - (self.m[3][0] * self.m[2][1] * self.m[1][3]) 
                    - (self.m[2][0] * self.m[1][1] * self.m[3][3])
                    - (self.m[1][0] * self.m[3][1] * self.m[2][3]),

                    - (self.m[1][0] * self.m[2][1] * self.m[3][2])
                    - (self.m[2][0] * self.m[3][1] * self.m[1][2])
                    - (self.m[3][0] * self.m[1][1] * self.m[2][2]) + 
                    (self.m[3][0] * self.m[2][1] * self.m[1][2]) + 
                    (self.m[2][0] * self.m[1][1] * self.m[3][2]) + 
                    (self.m[1][0] * self.m[3][1] * self.m[2][2])
                ], 
                [
                    - (self.m[0][1] * self.m[2][2] * self.m[3][3])
                    - (self.m[2][1] * self.m[3][2] * self.m[0][3])
                    - (self.m[3][1] * self.m[0][2] * self.m[2][3]) +
                    (self.m[3][1] * self.m[2][2] * self.m[0][3]) +
                    (self.m[2][1] * self.m[0][2] * self.m[3][3]) +
                    (self.m[0][1] * self.m[3][2] * self.m[2][3]),

                    (self.m[0][0] * self.m[2][2] * self.m[3][3]) +
                    (self.m[2][0] * self.m[3][2] * self.m[0][3]) +
                    (self.m[3][0] * self.m[0][2] * self.m[2][3])
                    - (self.m[3][0] * self.m[2][2] * self.m[0][3])
                    - (self.m[2][0] * self.m[0][2] * self.m[3][3])
                    - (self.m[0][0] * self.m[3][2] * self.m[2][3]),

                    - (self.m[0][0] * self.m[2][1] * self.m[3][3]) // ERROR
                    - (self.m[2][0] * self.m[3][1] * self.m[0][3])
                    - (self.m[3][0] * self.m[0][1] * self.m[2][3]) +
                    (self.m[3][0] * self.m[2][1] * self.m[0][3]) +
                    (self.m[2][0] * self.m[0][1] * self.m[3][3]) +
                    (self.m[0][0] * self.m[3][1] * self.m[2][3]), // FIX: self.m[1][3] to self.m[2][3] worked

                    (self.m[0][0] * self.m[2][1] * self.m[3][2]) +
                    (self.m[2][0] * self.m[3][1] * self.m[0][2]) +
                    (self.m[3][0] * self.m[0][1] * self.m[2][2])
                    - (self.m[3][0] * self.m[2][1] * self.m[0][2])
                    - (self.m[2][0] * self.m[0][1] * self.m[3][2])
                    - (self.m[0][0] * self.m[3][1] * self.m[2][2])
                ], 
                [
                    (self.m[0][1] * self.m[1][2] * self.m[3][3]) +
                    (self.m[1][1] * self.m[3][2] * self.m[0][3]) +
                    (self.m[3][1] * self.m[0][2] * self.m[1][3])
                    - (self.m[3][1] * self.m[1][2] * self.m[0][3])
                    - (self.m[1][1] * self.m[0][2] * self.m[3][3])
                    - (self.m[0][1] * self.m[3][2] * self.m[1][3]),

                    - (self.m[0][0] * self.m[1][2] * self.m[3][3])
                    - (self.m[1][0] * self.m[3][2] * self.m[0][3])
                    - (self.m[3][0] * self.m[0][2] * self.m[1][3]) +
                    (self.m[3][0] * self.m[1][2] * self.m[0][3]) +
                    (self.m[1][0] * self.m[0][2] * self.m[3][3]) +
                    (self.m[0][0] * self.m[3][2] * self.m[1][3]),

                    (self.m[0][0] * self.m[1][1] * self.m[3][3]) + 
                    (self.m[1][0] * self.m[3][1] * self.m[0][3]) + 
                    (self.m[3][0] * self.m[0][1] * self.m[1][3])
                    - (self.m[3][0] * self.m[1][1] * self.m[0][3]) 
                    - (self.m[1][0] * self.m[0][1] * self.m[3][3])
                    - (self.m[0][0] * self.m[3][1] * self.m[1][3]),

                    - (self.m[0][0] * self.m[1][1] * self.m[3][2])
                    - (self.m[1][0] * self.m[3][1] * self.m[0][2])
                    - (self.m[3][0] * self.m[0][1] * self.m[1][2]) + 
                    (self.m[3][0] * self.m[1][1] * self.m[0][2]) + 
                    (self.m[1][0] * self.m[0][1] * self.m[3][2]) + 
                    (self.m[0][0] * self.m[3][1] * self.m[1][2])
                ], 
                [
                    - (self.m[0][1] * self.m[1][2] * self.m[2][3])
                    - (self.m[1][1] * self.m[2][2] * self.m[0][3])
                    - (self.m[2][1] * self.m[0][2] * self.m[1][3]) +
                    (self.m[2][1] * self.m[1][2] * self.m[0][3]) +
                    (self.m[1][1] * self.m[0][2] * self.m[2][3]) +
                    (self.m[0][1] * self.m[2][2] * self.m[1][3]),

                    (self.m[0][0] * self.m[1][2] * self.m[2][3]) +
                    (self.m[1][0] * self.m[2][2] * self.m[0][3]) +
                    (self.m[2][0] * self.m[0][2] * self.m[1][3])
                    - (self.m[2][0] * self.m[1][2] * self.m[0][3])
                    - (self.m[1][0] * self.m[0][2] * self.m[2][3])
                    - (self.m[0][0] * self.m[2][2] * self.m[1][3]),

                    - (self.m[0][0] * self.m[1][1] * self.m[2][3])
                    - (self.m[1][0] * self.m[2][1] * self.m[0][3])
                    - (self.m[2][0] * self.m[0][1] * self.m[1][3]) +
                    (self.m[2][0] * self.m[1][1] * self.m[0][3]) +
                    (self.m[1][0] * self.m[0][1] * self.m[2][3]) +
                    (self.m[0][0] * self.m[2][1] * self.m[1][3]),

                    (self.m[0][0] * self.m[1][1] * self.m[2][2]) +
                    (self.m[1][0] * self.m[2][1] * self.m[0][2]) +
                    (self.m[2][0] * self.m[0][1] * self.m[1][2]) 
                    - (self.m[2][0] * self.m[1][1] * self.m[0][2])
                    - (self.m[1][0] * self.m[0][1] * self.m[2][2])
                    - (self.m[0][0] * self.m[2][1] * self.m[1][2])
                ]
            ],
        }
    }
    
    /// Get the adjugate matrix of a 3x3 matrix.
    pub fn adjugate(self) -> Self
    where T: Mul<T, Output = T> + Neg<Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Copy {
        // The adjugate matrix is equal to cofactor after it is transposed.
        self.cofactor().transpose()
    }

    /// Get the inverse of a 4x4 matrix.
    ///
    /// # Arguments
    /// * `self` - The matrix the function is being called for
    ///
    /// Inserts a value computed from `f` into the option if it is [`None`], then
    /// returns a mutable reference to the contained value.
    ///
    /// # Examples
    /// 
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let x = Matrix4::new(
    ///   [1.0,12.0,13.0,14.0],
    ///   [21.0,2.0,23.0,24.0],
    ///   [31.0,32.0,3.0,34.0],
    ///   [41.0,42.0,43.0,4.0]);
    /// let x_inverse = x.inverse();
    /// assert_eq!(x_inverse, Matrix4::new(
    ///   [-57.0/932.0,17.0/932.0,53.0/4660.0,37.0/4660.0],
    ///   [139.0/4660.0,-8.0/233.0,51.0/4660.0,2.0/233.0],
    ///   [25.0/932.0,69.0/4660.0,-21.0/932.0,41.0/4660.0],
    ///   [59.0/2330.0,67.0/4660.0,5.0/466.0,-15.0/932.0])); // TODO is this transposed should not be !!!!!!!!!!!!!!!
    /// ```
    ///
    pub fn inverse(&self) -> Matrix4<T>
    where T: Add<T, Output = T> + Mul<T, Output = T> + Sub<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Copy { // TODO result
        let det: T = self.determinant();

        // Adjugate matrix
        let adjugate_mat = self.adjugate();

        adjugate_mat / det
    }
}

impl<T> Construct<T> for Matrix4<T> where T: Construct<T> {}
impl<T> SqMatrix<T, Vector4<T>> for Matrix4<T> where T: Construct<T> {}

impl<T> fmt::Debug for Matrix4<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,"\n[ {:?}\t {:?}\t {:?}\t {:?}\t ]
                  \n[ {:?}\t {:?}\t {:?}\t {:?}\t ]
                  \n[ {:?}\t {:?}\t {:?}\t {:?}\t ]
                  \n[ {:?}\t {:?}\t {:?}\t {:?}\t ]",
                  self.m[0][0], self.m[0][1], self.m[0][2], self.m[0][3],
                  self.m[1][0], self.m[1][1], self.m[1][2], self.m[1][3],
                  self.m[2][0], self.m[2][1], self.m[2][2], self.m[2][3],
                  self.m[3][0], self.m[3][1], self.m[3][2], self.m[3][3])
    }
}

impl<T> Add for Matrix4<T>
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix4<T> + T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat0 = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let mat1 = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 2, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let expected = Matrix4::new(
    ///     [2, 4, 6, 8],
    ///     [6, 8, 3, 4],
    ///     [8, 2, 4, 6],
    ///     [4, 6, 8, 2]
    /// );
    /// assert_eq!(mat0 + mat1, expected);
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        Matrix4{
            m: [
                [self.m[0][0] + rhs.m[0][0], self.m[0][1] + rhs.m[0][1], self.m[0][2] + rhs.m[0][2], self.m[0][3] + rhs.m[0][3]],
                [self.m[1][0] + rhs.m[1][0], self.m[1][1] + rhs.m[1][1], self.m[1][2] + rhs.m[1][2], self.m[1][3] + rhs.m[1][3]],
                [self.m[2][0] + rhs.m[2][0], self.m[2][1] + rhs.m[2][1], self.m[2][2] + rhs.m[2][2], self.m[2][3] + rhs.m[2][3]],
                [self.m[3][0] + rhs.m[3][0], self.m[3][1] + rhs.m[3][1], self.m[3][2] + rhs.m[3][2], self.m[3][3] + rhs.m[3][3]],
            ]
        }
    }
}

impl<T> Add<T> for Matrix4<T>
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix4<T> + T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let val = 2;
    /// let expected = Matrix4::new(
    ///     [3, 4, 5, 6],
    ///     [5, 6, 3, 4],
    ///     [6, 3, 4, 5],
    ///     [4, 5, 6, 3]
    /// );
    /// assert_eq!(mat + val, expected);
    /// ```
    fn add(self, rhs: T) -> Self::Output {
        Matrix4 {
            m: [
                [self.m[0][0] + rhs, self.m[0][1] + rhs, self.m[0][2] + rhs, self.m[0][3] + rhs],
                [self.m[1][0] + rhs, self.m[1][1] + rhs, self.m[1][2] + rhs, self.m[1][3] + rhs],
                [self.m[2][0] + rhs, self.m[2][1] + rhs, self.m[2][2] + rhs, self.m[2][3] + rhs],
                [self.m[3][0] + rhs, self.m[3][1] + rhs, self.m[3][2] + rhs, self.m[3][3] + rhs],
            ]
        }
    }
}

impl<T> Sub for Matrix4<T>
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix4<T> - T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat0 = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let mat1 = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 2, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let expected = Matrix4::new(
    ///     [0, 0,  0, 0],
    ///     [0, 0, -1, 0],
    ///     [0, 0,  0, 0],
    ///     [0, 0,  0, 0]
    /// );
    /// assert_eq!(mat0 - mat1, expected);
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        Matrix4 {
            m: [
                [self.m[0][0] - rhs.m[0][0], self.m[0][1] - rhs.m[0][1], self.m[0][2] - rhs.m[0][2], self.m[0][3] - rhs.m[0][3]],
                [self.m[1][0] - rhs.m[1][0], self.m[1][1] - rhs.m[1][1], self.m[1][2] - rhs.m[1][2], self.m[1][3] - rhs.m[1][3]],
                [self.m[2][0] - rhs.m[2][0], self.m[2][1] - rhs.m[2][1], self.m[2][2] - rhs.m[2][2], self.m[2][3] - rhs.m[2][3]],
                [self.m[3][0] - rhs.m[3][0], self.m[3][1] - rhs.m[3][1], self.m[3][2] - rhs.m[3][2], self.m[3][3] - rhs.m[3][3]],
            ]
        }
    }
}

impl<T> Sub<T> for Matrix4<T>
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix4<T> - T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let val = 2;
    /// let expected = Matrix4::new(
    ///     [-1,  0,  1,  2],
    ///     [ 1,  2, -1,  0],
    ///     [ 2, -1,  0,  1],
    ///     [ 0,  1,  2, -1]
    /// );
    /// assert_eq!(mat - val, expected);
    /// ```
    fn sub(self, rhs: T) -> Self::Output {
        Matrix4 {
            m: [
                [self.m[0][0] - rhs, self.m[0][1] - rhs, self.m[0][2] - rhs, self.m[0][3] - rhs],
                [self.m[1][0] - rhs, self.m[1][1] - rhs, self.m[1][2] - rhs, self.m[1][3] - rhs],
                [self.m[2][0] - rhs, self.m[2][1] - rhs, self.m[2][2] - rhs, self.m[2][3] - rhs],
                [self.m[3][0] - rhs, self.m[3][1] - rhs, self.m[3][2] - rhs, self.m[3][3] - rhs],
            ]
        }
    }
}

impl<T> Mul for Matrix4<T> where T: Add<T, Output = T> + Mul<T, Output = T> + Copy{
    type Output = Self;

    /// Matrix4 * Matrix4 = Matrix4
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat0 = Matrix4::new(
    ///     [4, 6, 6, 0],
    ///     [6, 4, 5, 0],
    ///     [5, 6, 4, 0],
    ///     [0, 0, 0, 1] 
    /// );
    /// let mat1 = Matrix4::new(
    ///     [1, 2, 3, 0],
    ///     [2, 3, 1, 0],
    ///     [2, 1, 2, 0],
    ///     [0, 0, 0, 1] 
    /// );
    /// let expected = Matrix4::new(
    ///     [28, 32, 30, 0],
    ///     [24, 29, 32, 0],
    ///     [25, 32, 29, 0],
    ///     [0, 0, 0, 1] 
    /// );
    /// assert_eq!(mat0 * mat1, expected);
    /// ```
    fn mul(self, rhs: Self) -> Self::Output {
        Matrix4 {
            m: [
                [
                    self.m[0][0] * rhs.m[0][0]
                        + self.m[0][1] * rhs.m[1][0]
                        + self.m[0][2] * rhs.m[2][0]
                        + self.m[0][3] * rhs.m[3][0],
                    self.m[0][0] * rhs.m[0][1]
                        + self.m[0][1] * rhs.m[1][1]
                        + self.m[0][2] * rhs.m[2][1]
                        + self.m[0][3] * rhs.m[3][1],
                    self.m[0][0] * rhs.m[0][2]
                        + self.m[0][1] * rhs.m[1][2]
                        + self.m[0][2] * rhs.m[2][2]
                        + self.m[0][3] * rhs.m[3][2],
                    self.m[0][0] * rhs.m[0][3]
                        + self.m[0][1] * rhs.m[1][3]
                        + self.m[0][2] * rhs.m[2][3]
                        + self.m[0][3] * rhs.m[3][3],
                ],
                [
                    self.m[1][0] * rhs.m[0][0]
                        + self.m[1][1] * rhs.m[1][0]
                        + self.m[1][2] * rhs.m[2][0]
                        + self.m[1][3] * rhs.m[3][0],
                    self.m[1][0] * rhs.m[0][1]
                        + self.m[1][1] * rhs.m[1][1]
                        + self.m[1][2] * rhs.m[2][1]
                        + self.m[1][3] * rhs.m[3][1],
                    self.m[1][0] * rhs.m[0][2]
                        + self.m[1][1] * rhs.m[1][2]
                        + self.m[1][2] * rhs.m[2][2]
                        + self.m[1][3] * rhs.m[2][3],
                    self.m[1][0] * rhs.m[0][3]
                        + self.m[1][1] * rhs.m[1][3]
                        + self.m[1][2] * rhs.m[2][3]
                        + self.m[1][3] * rhs.m[3][3],
                ],
                [
                    self.m[2][0] * rhs.m[0][0]
                        + self.m[2][1] * rhs.m[1][0]
                        + self.m[2][2] * rhs.m[2][0]
                        + self.m[2][3] * rhs.m[3][0],
                    self.m[2][0] * rhs.m[0][1]
                        + self.m[2][1] * rhs.m[1][1]
                        + self.m[2][2] * rhs.m[2][1]
                        + self.m[2][3] * rhs.m[3][1],
                    self.m[2][0] * rhs.m[0][2]
                        + self.m[2][1] * rhs.m[1][2]
                        + self.m[2][2] * rhs.m[2][2]
                        + self.m[2][3] * rhs.m[3][2],
                    self.m[2][0] * rhs.m[0][3]
                        + self.m[2][1] * rhs.m[1][3]
                        + self.m[2][2] * rhs.m[2][3]
                        + self.m[2][3] * rhs.m[3][3],
                ],
                [
                    self.m[3][0] * rhs.m[0][0]
                        + self.m[3][1] * rhs.m[1][0]
                        + self.m[3][2] * rhs.m[2][0]
                        + self.m[3][3] * rhs.m[3][0],
                    self.m[3][0] * rhs.m[0][1]
                        + self.m[3][1] * rhs.m[1][1]
                        + self.m[3][2] * rhs.m[2][1]
                        + self.m[3][3] * rhs.m[3][1],
                    self.m[3][0] * rhs.m[0][2]
                        + self.m[3][1] * rhs.m[1][2]
                        + self.m[3][2] * rhs.m[2][2]
                        + self.m[3][3] * rhs.m[3][2],
                    self.m[3][0] * rhs.m[0][3]
                        + self.m[3][1] * rhs.m[1][3]
                        + self.m[3][2] * rhs.m[2][3]
                        + self.m[3][3] * rhs.m[3][3],
                ],
            ],
        }
    }
}

impl<T> Mul<Vector4<T>> for Matrix4<T> where T: Add<T, Output = T> + Mul<T, Output = T> + Copy{
    type Output = Vector4<T>;

    /// Matrix4 * Vec3d = type Vec3d
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector4::Vector4;
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let vec = Vector4::new(1, 4, 2, 3);
    /// assert_eq!(mat * vec, Vector4::new(27, 27, 21, 25));
    /// ```
    fn mul(self, rhs: Vector4<T>) -> Self::Output {
        Vector4::new(
            rhs.0 * self.m[0][0] + rhs.1 * self.m[0][1] + rhs.2 * self.m[0][2] + rhs.3 * self.m[0][3],
            rhs.0 * self.m[1][0] + rhs.1 * self.m[1][1] + rhs.2 * self.m[1][2] + rhs.3 * self.m[1][3],
            rhs.0 * self.m[2][0] + rhs.1 * self.m[2][1] + rhs.2 * self.m[2][2] + rhs.3 * self.m[2][3],
            rhs.0 * self.m[3][0] + rhs.1 * self.m[3][1] + rhs.2 * self.m[3][2] + rhs.3 * self.m[3][3],
        )
    }
}

impl<T> Mul<T> for Matrix4<T> where T: Mul<T, Output = T> + Copy{
    type Output = Matrix4<T>;

    /// Matrix4<T> * T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let val = 2;
    /// let expected = Matrix4::new(
    ///     [2, 4, 6, 8],
    ///     [6, 8, 2, 4],
    ///     [8, 2, 4, 6],
    ///     [4, 6, 8, 2]
    /// );
    /// assert_eq!(mat * val, expected);
    /// ```
    fn mul(self, rhs: T) -> Self::Output {
        Matrix4 {
            m: [
                [self.m[0][0] * rhs, self.m[0][1] * rhs, self.m[0][2] * rhs, self.m[0][3] * rhs],
                [self.m[1][0] * rhs, self.m[1][1] * rhs, self.m[1][2] * rhs, self.m[1][3] * rhs],
                [self.m[2][0] * rhs, self.m[2][1] * rhs, self.m[2][2] * rhs, self.m[2][3] * rhs],
                [self.m[3][0] * rhs, self.m[3][1] * rhs, self.m[3][2] * rhs, self.m[3][3] * rhs],
            ],
        }
    }
}

impl<T> Div for Matrix4<T>
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Copy {
    type Output = Self;

    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector4::Vector4;
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat0 = Matrix4::new(
    ///     [4, 5, 6, 0],
    ///     [6, 4, 5, 0],
    ///     [5, 6, 4, 0],
    ///     [0, 0, 0, 1],
    /// );
    /// let mat1 = Matrix4::new(
    ///     [1, 2, 3, 0],
    ///     [2, 3, 1, 0],
    ///     [3, 1, 2, 0],
    ///     [0, 0, 0, 1]
    /// );
    /// let expected = Matrix4::new(
    ///     [1, 0, 0, 0],
    ///     [0, 0, 1, 0],
    ///     [0, 1, 0, 0],
    ///     [0, 0, 0, 1],
    /// );
    /// assert_eq!(mat0 / mat1, expected);
    /// ```
    fn div(self, rhs: Self) -> Self::Output {
        // self * rhs.inverse()

        let determinant = rhs.determinant();

        let adjugate_mat = rhs.adjugate();

        // first multiply to mitigate error with non Floatingpoint types
        self * adjugate_mat / determinant
    }
}

impl<T> Div<T> for Matrix4<T> where T: Div<T, Output = T> + Copy {
    type Output = Matrix4<T>;

    /// Matrix4<T> / T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let val = 2;
    /// let expected = Matrix4::new(
    ///     [0, 1, 1, 2],
    ///     [1, 2, 0, 1],
    ///     [2, 0, 1, 1],
    ///     [1, 1, 2, 0]
    /// );
    /// assert_eq!(mat / val, expected);
    /// ```
    fn div(self, rhs: T) -> Self::Output {
        Matrix4 {
            m: [
                [self.m[0][0] / rhs, self.m[0][1] / rhs, self.m[0][2] / rhs, self.m[0][3] / rhs],
                [self.m[1][0] / rhs, self.m[1][1] / rhs, self.m[1][2] / rhs, self.m[1][3] / rhs],
                [self.m[2][0] / rhs, self.m[2][1] / rhs, self.m[2][2] / rhs, self.m[2][3] / rhs],
                [self.m[3][0] / rhs, self.m[3][1] / rhs, self.m[3][2] / rhs, self.m[3][3] / rhs],
            ],
        }
    }
}


impl<T> Rem for Matrix4<T>
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix4<T> % Matrix4<T>
    /// = Matrix4<T>
    /// Finds the remainder after an ELEMENT WISE division
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat0 = Matrix4::new(
    ///     [4, 5, 6, 3],
    ///     [6, 4, 5, 2],
    ///     [5, 6, 4, 1],
    ///     [8, 1, 3, 2],
    /// );
    /// let mat1 = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [2, 3, 1, 5],
    ///     [3, 1, 2, 6],
    ///     [3, 9, 8, 7],
    /// );
    /// let expected = Matrix4::new(
    ///     [0, 1, 0, 3],
    ///     [0, 1, 0, 2],
    ///     [2, 0, 0, 1],
    ///     [2, 1, 3, 2],
    /// );
    /// assert_eq!(mat0 % mat1, expected);
    /// ```
    fn rem(self, rhs: Self) -> Self::Output {
        Matrix4::new(
            [self.m[0][0] % rhs.m[0][0], self.m[0][1] % rhs.m[0][1], self.m[0][2] % rhs.m[0][2], self.m[0][3] % rhs.m[0][3]],
            [self.m[1][0] % rhs.m[1][0], self.m[1][1] % rhs.m[1][1], self.m[1][2] % rhs.m[1][2], self.m[1][3] % rhs.m[1][3]],
            [self.m[2][0] % rhs.m[2][0], self.m[2][1] % rhs.m[2][1], self.m[2][2] % rhs.m[2][2], self.m[2][3] % rhs.m[2][3]],
            [self.m[3][0] % rhs.m[3][0], self.m[3][1] % rhs.m[3][1], self.m[3][2] % rhs.m[3][2], self.m[3][3] % rhs.m[3][3]],
        )
    }
}

impl<T> Rem<T> for Matrix4<T> where T: Rem<T, Output = T> + Copy {
    type Output = Matrix4<T>;

    /// Matrix4<T> % T
    /// = Matrix4<T>
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let val = 2;
    /// let expected = Matrix4::new(
    ///     [1, 0, 1, 0],
    ///     [1, 0, 1, 0],
    ///     [0, 1, 0, 1],
    ///     [0, 1, 0, 1]
    /// );
    /// assert_eq!(mat % val, expected);
    /// ```
    /// ```should_panic
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// mat % 0;
    /// ```
    fn rem(self, rhs: T) -> Self::Output {
        Matrix4 {
            m: [
                [self.m[0][0] % rhs, self.m[0][1] % rhs, self.m[0][2] % rhs, self.m[0][3] % rhs],
                [self.m[1][0] % rhs, self.m[1][1] % rhs, self.m[1][2] % rhs, self.m[1][3] % rhs],
                [self.m[2][0] % rhs, self.m[2][1] % rhs, self.m[2][2] % rhs, self.m[2][3] % rhs],
                [self.m[3][0] % rhs, self.m[3][1] % rhs, self.m[3][2] % rhs, self.m[3][3] % rhs],
            ],
        }
    }
}


impl<T> Neg for Matrix4<T>
where T: Neg<Output = T> + Copy {
    type Output = Self;
    
    /// -Matrix4<T>
    /// = Matrix4<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::matrix4::Matrix4;
    /// let mat = Matrix4::new(
    ///     [1, 2, 3, 4],
    ///     [3, 4, 1, 2],
    ///     [4, 1, 2, 3],
    ///     [2, 3, 4, 1]
    /// );
    /// let expected = Matrix4::new(
    ///     [-1, -2, -3, -4],
    ///     [-3, -4, -1, -2],
    ///     [-4, -1, -2, -3],
    ///     [-2, -3, -4, -1]
    /// );
    /// assert_eq!(-mat, expected);
    /// ```
    fn neg(self) -> Self::Output {
        Matrix4{
            m: [
                [-self.m[0][0], -self.m[0][1], -self.m[0][2], -self.m[0][3]],
                [-self.m[1][0], -self.m[1][1], -self.m[1][2], -self.m[1][3]],
                [-self.m[2][0], -self.m[2][1], -self.m[2][2], -self.m[2][3]],
                [-self.m[3][0], -self.m[3][1], -self.m[3][2], -self.m[3][3]],
            ],
        }
    }
}

impl<T> Zero for Matrix4<T> 
where T: Zero + Copy {
    const ZERO: Self = Matrix4 {
        m: [
            [T::ZERO, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ZERO]
        ]
    };
}

impl<T> One for Matrix4<T> 
where T: Zero + One {
    /// The identity matrix
    const ONE: Self = Matrix4 {
        m: [
            [T::ONE , T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ONE , T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ONE , T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE ]
        ]
    };
}

impl<T> Two for Matrix4<T> 
where T: Zero + Two {
    /// The identity matrix
    const TWO: Self = Matrix4 {
        m: [
            [T::TWO , T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::TWO , T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::TWO , T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::TWO ]
        ]
    };
}


impl<T> From<Matrix4<T>> for [[T; 4]; 4] {
    fn from(other: Matrix4<T>) -> [[T; 4]; 4] {
        other.m
    }
}

impl<T> F32Fmt for Matrix4<T> 
where T: F32Fmt + Copy { 
    type F32Fmt = Matrix4<T::F32Fmt>;
    #[inline] 
    fn intoF32Fmt(self) -> Self::F32Fmt {
        let m = &self.m;
        Matrix4{
            m: [
                [m[0][0].intoF32Fmt(), m[0][1].intoF32Fmt(), m[0][2].intoF32Fmt(), m[0][3].intoF32Fmt()],
                [m[1][0].intoF32Fmt(), m[1][1].intoF32Fmt(), m[1][2].intoF32Fmt(), m[1][3].intoF32Fmt()],
                [m[2][0].intoF32Fmt(), m[2][1].intoF32Fmt(), m[2][2].intoF32Fmt(), m[2][3].intoF32Fmt()],
                [m[3][0].intoF32Fmt(), m[3][1].intoF32Fmt(), m[3][2].intoF32Fmt(), m[3][3].intoF32Fmt()],
            ]
        }
    } 
    #[inline]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { 
        let m = &f32_fmt.m;
        Matrix4{
            m: [
                [T::fromF32Fmt(m[0][0]), T::fromF32Fmt(m[0][1]), T::fromF32Fmt(m[0][2]), T::fromF32Fmt(m[0][2])],
                [T::fromF32Fmt(m[1][0]), T::fromF32Fmt(m[1][1]), T::fromF32Fmt(m[1][2]), T::fromF32Fmt(m[1][2])],
                [T::fromF32Fmt(m[2][0]), T::fromF32Fmt(m[2][1]), T::fromF32Fmt(m[2][2]), T::fromF32Fmt(m[2][2])],
                [T::fromF32Fmt(m[3][0]), T::fromF32Fmt(m[3][1]), T::fromF32Fmt(m[3][2]), T::fromF32Fmt(m[3][2])],
            ]
        }
    }

    fn sqrt(self) -> Self {
        todo!()
    }

    fn cbrt(self) -> Self { 
        todo!() 
    }

    fn f32_const_mul(self, constant: f32) -> Self {
        Matrix4{
            m: [
                [self.m[0][0].f32_const_mul(constant), self.m[0][1].f32_const_mul(constant), self.m[0][2].f32_const_mul(constant), self.m[0][3].f32_const_mul(constant)],
                [self.m[1][0].f32_const_mul(constant), self.m[1][1].f32_const_mul(constant), self.m[1][2].f32_const_mul(constant), self.m[1][3].f32_const_mul(constant)],
                [self.m[2][0].f32_const_mul(constant), self.m[2][1].f32_const_mul(constant), self.m[2][2].f32_const_mul(constant), self.m[2][3].f32_const_mul(constant)],
                [self.m[3][0].f32_const_mul(constant), self.m[3][1].f32_const_mul(constant), self.m[3][2].f32_const_mul(constant), self.m[3][3].f32_const_mul(constant)],
            ]
        }
    }

    fn sin_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn cos_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn tan_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn asin_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn acos_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn atan_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn atan2_mul(self, _other: Self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn sinh_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn cosh_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }

    fn tanh_mul(self, _mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        todo!()
    }
}

impl<T> SignOps for Matrix4<T> {
    fn ptcopysign(self, _sign: Self) -> Self {
        todo!()
    }

    fn ptsignum(self) -> i8 {
        todo!()
    }

    fn abs(self) -> Self {
        todo!()
    }
}