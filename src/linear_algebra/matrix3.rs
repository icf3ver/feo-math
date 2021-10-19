use crate::{F32Fmt, Primitive, SignOps, /*Sum,*/ Three, Two, polynomial::CubicEquation};

use {
    super::SqMatrix,
    crate::{
        Construct, 
        multivariate::{
            two_unknown::Equa2,
            with_var::WithVar
        },
        linear_algebra::{
            matrix4::Matrix4, 
            vector3::Vector3
        },
        One, Zero
    },
    std::{
        fmt, 
        ops::{Add, Div, Rem, Mul, Neg, Sub}
    }
};

/// A 3 x 3 Matrix with row major order
#[derive(PartialEq, Clone, Copy)]
pub struct Matrix3<T>{pub m:[[T; 3]; 3]}

impl<T> Matrix3<T> {
    pub fn new(
            row0: [T; 3], 
            row1: [T; 3], 
            row2: [T; 3]) -> Self {
        Matrix3 {
            m: [
                row0,
                row1,
                row2
            ]
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_val (
            r0c0: T, r0c1: T, r0c2: T,
            r1c0: T, r1c1: T, r1c2: T,
            r2c0: T, r2c1: T, r2c2: T) -> Self {
        Matrix3 {
            m: [
                [r0c0, r0c1, r0c2],
                [r1c0, r1c1, r1c2],
                [r2c0, r2c1, r2c2]
            ]
        }
    }

    /// return the identity
    pub fn identity() -> Self
    where T: Zero + One {
        Self::ONE
    }

    /// the trace ...
    pub fn trace(&self) -> T
    where T: Add<T, Output = T> + Copy {
        self.m[0][0] + self.m[1][1] + self.m[2][2]
    }

    /// Calculate the determinant
    ///
    /// # Arguments
    ///
    /// * `self` - The Matrix3 the function was called for
    ///
    #[allow(dead_code)]
    pub fn determinant(&self) -> T 
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Copy {
        // matrix (a, b, c)
        //        (d, e, f)
        //        (g, h, i)
        // determinant of a 3x3 matrix is |A| = a(ei − fh)
        //                                    − b(di − fg)
        //                                    + c(dh − eg)

        self.m[0][0] * ((self.m[1][1] * self.m[2][2]) - (self.m[1][2] * self.m[2][1]))
            - self.m[0][1] * ((self.m[1][0] * self.m[2][2]) - (self.m[1][2] * self.m[2][0]))
            + self.m[0][2] * ((self.m[1][0] * self.m[2][1]) - (self.m[1][1] * self.m[2][0]))
    }

    /// Return a transposed matrix
    /// transpose
    /// matrix (a, b, c)
    ///        (d, e, f)
    ///        (g, h, i)
    ///
    /// matrix (a, d, g)
    ///        (b, e, h)
    ///        (c, f, i)
    /// ```rust
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let transposed = Matrix3::new(
    ///     [0.0,2.0,2.0],
    ///     [1.0,0.0,2.0],
    ///     [1.0,1.0,0.0]
    /// ).transpose();
    /// let result = Matrix3::new(
    ///     [0.0,1.0,1.0],
    ///     [2.0,0.0,1.0],
    ///     [2.0,2.0,0.0]
    /// );
    /// assert_eq!(transposed, result);
    /// ```
    pub fn transpose(&self) -> Matrix3<T>
    where T: Copy {
        Matrix3 {
            m: [
                [self.m[0][0], self.m[1][0], self.m[2][0]],
                [self.m[0][1], self.m[1][1], self.m[2][1]],
                [self.m[0][2], self.m[1][2], self.m[2][2]],
            ],
        }
    }

    /// Rule:
    /// Mv = λv
    /// where M = a matrix
    ///       λ = an eigenvalue
    ///       v = a eigenvector
    ///
    /// as a 3x3 Matrix can have  
    ///
    /// if v != <0, 0, 0> then:
    ///     ... // self expl todo
    ///     |M - Iλ| = 0
    ///     ... // in fn
    pub fn eigen_value(&self) -> [Option<T>; 3] 
    where T: F32Fmt<F32Fmt = f32> + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + SignOps + F32Fmt + Zero + One + Three + PartialEq + Copy
        + Rem<T, Output = T> + Two + Primitive + fmt::Debug {
        let m = self.m;
        // 0 = (m[0][0] - λ)((m[1][1] - λ)(m[2][2] - λ) - m[1][2]m[2][1])
        //     - m[0][1](m[1][0](m[2][2] - λ) - m[1][2]m[2][0])
        //     + m[0][2](m[1][0]m[2][1] - (m[1][1] - λ)m[2][0])

        // 0 = (m[0][0] - λ)(m[1][1] - λ)(m[2][2] - λ) - (m[0][0] - λ)m[1][2]m[2][1]
        //     - m[0][1]m[1][0](m[2][2] - λ) + m[0][1]m[1][2]m[2][0]
        //     + m[0][2]m[1][0]m[2][1] - m[0][2](m[1][1] - λ)m[2][0]

        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] = 
        //      (m[0][0] - λ)(m[1][1] - λ)(m[2][2] - λ) 
        //      - (m[0][0] - λ)m[1][2]m[2][1]
        //      - m[0][1]m[1][0](m[2][2] - λ)
        //      - m[0][2](m[1][1] - λ)m[2][0]

        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] = 
        //      (m[0][0] - λ)(m[1][1] - λ)(m[2][2] - λ) 
        //      - m[0][0]m[1][2]m[2][1] + m[1][2]m[2][1]λ
        //      - m[0][1]m[1][0]m[2][2] + m[0][1]m[1][0]λ
        //      - m[0][2]m[1][1]m[2][0] + m[0][2]m[2][0]λ

        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] + m[0][0]m[1][2]m[2][1] + m[0][1]m[1][0]m[2][2] + m[0][2]m[1][1]m[2][0] = 
        //      (m[0][0] - λ)(m[1][1] - λ)(m[2][2] - λ) 
        //      + m[1][2]m[2][1]λ
        //      + m[0][1]m[1][0]λ
        //      + m[0][2]m[2][0]λ

        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] + m[0][0]m[1][2]m[2][1] + m[0][1]m[1][0]m[2][2] + m[0][2]m[1][1]m[2][0] = 
        //      (m[0][0]m[1][1] - m[0][0]λ - m[1][1]λ + λ^2)(m[2][2] - λ) 
        //      + m[1][2]m[2][1]λ
        //      + m[0][1]m[1][0]λ
        //      + m[0][2]m[2][0]λ
        
        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] + m[0][0]m[1][2]m[2][1] + m[0][1]m[1][0]m[2][2] + m[0][2]m[1][1]m[2][0] = 
        //      (m[0][0]m[1][1]m[2][2] - m[0][0]m[2][2]λ - m[1][1]m[2][2]λ + m[2][2]λ^2 - m[0][0]m[1][1]λ + m[0][0]λ^2 + m[1][1]λ^2 - λ^3)
        //      + m[1][2]m[2][1]λ
        //      + m[0][1]m[1][0]λ
        //      + m[0][2]m[2][0]λ
        
        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] + m[0][0]m[1][2]m[2][1] + m[0][1]m[1][0]m[2][2] + m[0][2]m[1][1]m[2][0] = 
        //      (m[0][0]m[1][1]m[2][2] - (m[0][0]m[2][2] + m[1][1]m[2][2] + m[0][0]m[1][1])λ + (m[2][2] + m[0][0] + m[1][1])λ^2 - λ^3)
        //      + m[1][2]m[2][1]λ
        //      + m[0][1]m[1][0]λ
        //      + m[0][2]m[2][0]λ
        
        // -m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] + m[0][0]m[1][2]m[2][1] + m[0][1]m[1][0]m[2][2] + m[0][2]m[1][1]m[2][0] - m[0][0]m[1][1]m[2][2] = 
        //      - (m[0][0]m[2][2] + m[1][1]m[2][2] + m[0][0]m[1][1])λ + (m[2][2] + m[0][0] + m[1][1])λ^2 - λ^3
        //      + m[1][2]m[2][1]λ
        //      + m[0][1]m[1][0]λ
        //      + m[0][2]m[2][0]λ

        // 0 =  - λ^3
        //      + (m[2][2] + m[0][0] + m[1][1])λ^2
        //      - (m[0][0]m[2][2] + m[1][1]m[2][2] + m[0][0]m[1][1] - m[1][2]m[2][1] - m[0][1]m[1][0] - m[0][2]m[2][0])λ
        //      - (-m[0][1]m[1][2]m[2][0] - m[0][2]m[1][0]m[2][1] + m[0][0]m[1][2]m[2][1] + m[0][1]m[1][0]m[2][2] + m[0][2]m[1][1]m[2][0] - m[0][0]m[1][1]m[2][2])
        
        // 0 =  - λ^3
        //      + (m[2][2] + m[0][0] + m[1][1])λ^2
        //      + (-m[0][0]m[2][2] - m[1][1]m[2][2] - m[0][0]m[1][1] + m[1][2]m[2][1] + m[0][1]m[1][0] + m[0][2]m[2][0])λ
        //      + (m[0][1]m[1][2]m[2][0] + m[0][2]m[1][0]m[2][1] - m[0][0]m[1][2]m[2][1] - m[0][1]m[1][0]m[2][2] - m[0][2]m[1][1]m[2][0] + m[0][0]m[1][1]m[2][2])
        
        // 0 =  - λ^3
        //      + (m[2][2] + m[0][0] + m[1][1])λ^2
        //      + (-m[0][0]m[2][2] - m[1][1]m[2][2] - m[0][0]m[1][1] + m[1][2]m[2][1] + m[0][1]m[1][0] + m[0][2]m[2][0])λ
        //      + (m[0][1]m[1][2]m[2][0] + m[0][2]m[1][0]m[2][1] - m[0][0]m[1][2]m[2][1] - m[0][1]m[1][0]m[2][2] - m[0][2]m[1][1]m[2][0] + m[0][0]m[1][1]m[2][2])

        let eigenvalues = CubicEquation::new(
            -T::ONE,
            m[2][2] + m[0][0] + m[1][1],
            -m[0][0] * m[2][2] - m[1][1] * m[2][2] - m[0][0] * m[1][1] + m[1][2] * m[2][1] + m[0][1] * m[1][0] + m[0][2] * m[2][0],
            m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1] - m[0][0] * m[1][2] * m[2][1] - m[0][1] * m[1][0] * m[2][2] - m[0][2] * m[1][1] * m[2][0] + m[0][0] * m[1][1] * m[2][2],
        ).solve(); // TOFIX
        
        let mut real_eigenvalues = eigenvalues.iter().map(|answer| match *answer {
            crate::imaginary::PossibleImag::Complex(r, _) => Some(r), // very incorrect
            crate::imaginary::PossibleImag::Imag(_) => None,
            crate::imaginary::PossibleImag::Real(r) => Some(r),
        });

        [real_eigenvalues.next().unwrap(), real_eigenvalues.next().unwrap(), real_eigenvalues.next().unwrap()]
    }

    pub fn eigen_vector(&self, eigen_value: T) -> Vector3<T>
    where T:  Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + fmt::Debug + Zero + One + Copy {
        let m = (*self - (Matrix3::<T>::identity() * eigen_value)).m;
        let equation_1: Equa2<T> = Equa2::new(m[0][0], WithVar{var: 'y', val: m[1][0]}, WithVar{var: 'z', val: m[2][0]});
        let equation_2: Equa2<T> = Equa2::new(m[0][1], WithVar{var: 'y', val: m[1][1]}, WithVar{var: 'z', val: m[2][1]});
        let y_z = equation_1.set_eq(equation_2); // rmv self make Equa2::set_eq()
        Vector3(T::ONE, y_z.0.1, y_z.1.1)
    }

    /// Get the cofactor matrix of a 3x3 matrix.
    #[allow(clippy::many_single_char_names)]
    pub fn cofactor(self) -> Self 
    where T: Mul<T, Output = T> + Neg<Output = T> + Sub<T, Output = T> + Copy {
        // Cofactor of a matrix
        // matrix (a, b, c)
        //        (d, e, f)
        //        (g, h, i)
        //
        // matrix ( |e, f|,  |d, f|,  |d, e|)
        //        ( |h, i|  -|g, i|   |g, h|)
        //
        //        ( |b, c|,  |a, c|,  |a, b|)
        //        (-|h, i|   |g, i|  -|g, h|)
        //
        //        ( |b, c|,  |a, c|,  |a, b|)
        //        ( |e, f|  -|d, f|   |d, e|)
        //
        // matrix (a, b)
        //        (c, d)
        // determinant of 2x2 matrix if |A| = (a * d) - (b * c)
        //
        // so
        //
        // matrix ( ((e * i) - (f * h)), -((d * i) - (f * g)),  ((d * h) - (e * g)))
        //        (-((b * i) - (c * h)),  ((a * i) - (c * g)), -((a * h) - (b * g)))
        //        ( ((b * f) - (c * e)), -((a * f) - (c * d)),  ((a * e) - (b * d)))

        let Matrix3{
            m: [
                [a, b, c],
                [d, e, f],
                [g, h, i]
            ]
        } = self;

        let a0_0 = (e * i) - (f * h);
        let b0_1 = -((d * i) - (f * g));
        let c0_2 = (d * h) - (e * g);
        let d1_0 = -((b * i) - (c * h));
        let e1_1 = (a * i) - (c * g);
        let f1_2 = -((a * h) - (b * g));
        let g2_0 = (b * f) - (c * e);
        let h2_1 = -((a * f) - (c * d));
        let i2_2 = (a * e) - (b * d);
        
        Matrix3 {
            m: [
                [a0_0, b0_1, c0_2], 
                [d1_0, e1_1, f1_2], 
                [g2_0, h2_1, i2_2]
            ],
        }
    }

    /// Get the adjugate matrix of a 3x3 matrix.
    pub fn adjugate(self) -> Self
    where T: Mul<T, Output = T> + Neg<Output = T> + Sub<T, Output = T> + Copy {
        // The adjugate matrix is equal to the transposed cofactor.
        self.cofactor().transpose()
    }

    /// Get the inverse of a 3x3 matrix.
    /// ```rust
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat0 = Matrix3::new(
    ///     [1, 1, 0],
    ///     [0, 1, 0],
    ///     [1, 0, 1],
    /// );
    /// let expected = Matrix3::new(
    ///     [ 1, -1,  0],
    ///     [ 0,  1,  0],
    ///     [-1,  1,  1],
    /// );
    /// assert_eq!(mat0.inverse(), expected);
    /// ```
    /// # Arguments
    /// * `self` - The matrix the function is being called for
    pub fn inverse(&self) -> Self 
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Copy {

        let determinant: T = self.determinant();

        let adjugate_mat = self.adjugate();

        adjugate_mat / determinant
    }
}

impl<T> Construct<T> for Matrix3<T> where T: Construct<T> {}
impl<T> SqMatrix<T, Vector3<T>> for Matrix3<T> where T: Construct<T> {}

impl<T> fmt::Debug for Matrix3<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,"\n[ {:?}\t {:?}\t {:?}\t ]
                  \n[ {:?}\t {:?}\t {:?}\t ]
                  \n[ {:?}\t {:?}\t {:?}\t ]",
                  self.m[0][0], self.m[0][1], self.m[0][2],
                  self.m[1][0], self.m[1][1], self.m[1][2],
                  self.m[2][0], self.m[2][1], self.m[2][2])
    }
}

impl<T> Add for Matrix3<T>
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> + Matrix3<T>
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat0 = Matrix3::new(
    ///     [4, 5, 6],
    ///     [6, 4, 5],
    ///     [5, 6, 4]
    /// );
    /// let mat1 = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let expected = Matrix3::new(
    ///     [5, 7, 9],
    ///     [8, 7, 6],
    ///     [8, 7, 6] 
    /// );
    /// assert_eq!(mat0 + mat1, expected);
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        Matrix3{
            m: [
                [self.m[0][0] + rhs.m[0][0], self.m[0][1] + rhs.m[0][1], self.m[0][2] + rhs.m[0][2]],
                [self.m[1][0] + rhs.m[1][0], self.m[1][1] + rhs.m[1][1], self.m[1][2] + rhs.m[1][2]],
                [self.m[2][0] + rhs.m[2][0], self.m[2][1] + rhs.m[2][1], self.m[2][2] + rhs.m[2][2]],
            ]
        }
    }
}

impl<T> Add<T> for Matrix3<T>
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> + T
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let val = 2;
    /// let expected = Matrix3::new(
    ///     [3, 4, 5],
    ///     [4, 5, 3],
    ///     [5, 3, 4]
    /// );
    /// assert_eq!(mat + val, expected);
    /// ```
    fn add(self, rhs: T) -> Self::Output {
        Matrix3{
            m: [
                [self.m[0][0] + rhs, self.m[0][1] + rhs, self.m[0][2] + rhs],
                [self.m[1][0] + rhs, self.m[1][1] + rhs, self.m[1][2] + rhs],
                [self.m[2][0] + rhs, self.m[2][1] + rhs, self.m[2][2] + rhs],
            ]
        }
    }
}

impl<T> Sub for Matrix3<T>
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> - Matrix3<T>
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat0 = Matrix3::new(
    ///     [4, 5, 6],
    ///     [6, 4, 5],
    ///     [5, 6, 4]
    /// );
    /// let mat1 = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let expected = Matrix3::new(
    ///     [3, 3, 3],
    ///     [4, 1, 4],
    ///     [2, 5, 2] 
    /// );
    /// assert_eq!(mat0 - mat1, expected);
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        Matrix3{
            m: [
                [self.m[0][0] - rhs.m[0][0], self.m[0][1] - rhs.m[0][1], self.m[0][2] - rhs.m[0][2]],
                [self.m[1][0] - rhs.m[1][0], self.m[1][1] - rhs.m[1][1], self.m[1][2] - rhs.m[1][2]],
                [self.m[2][0] - rhs.m[2][0], self.m[2][1] - rhs.m[2][1], self.m[2][2] - rhs.m[2][2]],
            ]
        }
    }
}

impl<T> Sub<T> for Matrix3<T>
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> - T
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let val = 2;
    /// let expected = Matrix3::new(
    ///     [-1,  0,  1],
    ///     [ 0,  1, -1],
    ///     [ 1, -1,  0]
    /// );
    /// assert_eq!(mat - val, expected);
    /// ```
    fn sub(self, rhs: T) -> Self::Output {
        Matrix3{
            m: [
                [self.m[0][0] - rhs, self.m[0][1] - rhs, self.m[0][2] - rhs],
                [self.m[1][0] - rhs, self.m[1][1] - rhs, self.m[1][2] - rhs],
                [self.m[2][0] - rhs, self.m[2][1] - rhs, self.m[2][2] - rhs],
            ]
        }
    }
}

impl<T> Mul for Matrix3<T>
where T: Add<T, Output = T> + Mul<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> * Matrix3<T>
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat0 = Matrix3::new(
    ///     [4, 5, 6],
    ///     [6, 4, 5],
    ///     [5, 6, 4]
    /// );
    /// let mat1 = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [2, 1, 2]
    /// );
    /// let expected = Matrix3::new(
    ///     [26, 29, 29],
    ///     [24, 29, 32],
    ///     [25, 32, 29] 
    /// );
    /// assert_eq!(mat0 * mat1, expected);
    /// ```
    fn mul(self, rhs: Self) -> Self::Output {
        Matrix3 {
            m: [
                [
                    self.m[0][0] * rhs.m[0][0]
                        + self.m[0][1] * rhs.m[1][0]
                        + self.m[0][2] * rhs.m[2][0],
                    self.m[0][0] * rhs.m[0][1]
                        + self.m[0][1] * rhs.m[1][1]
                        + self.m[0][2] * rhs.m[2][1],
                    self.m[0][0] * rhs.m[0][2]
                        + self.m[0][1] * rhs.m[1][2]
                        + self.m[0][2] * rhs.m[2][2],
                ],
                [
                    self.m[1][0] * rhs.m[0][0]
                        + self.m[1][1] * rhs.m[1][0]
                        + self.m[1][2] * rhs.m[2][0],
                    self.m[1][0] * rhs.m[0][1]
                        + self.m[1][1] * rhs.m[1][1]
                        + self.m[1][2] * rhs.m[2][1],
                    self.m[1][0] * rhs.m[0][2]
                        + self.m[1][1] * rhs.m[1][2]
                        + self.m[1][2] * rhs.m[2][2],
                ],
                [
                    self.m[2][0] * rhs.m[0][0]
                        + self.m[2][1] * rhs.m[1][0]
                        + self.m[2][2] * rhs.m[2][0],
                    self.m[2][0] * rhs.m[0][1]
                        + self.m[2][1] * rhs.m[1][1]
                        + self.m[2][2] * rhs.m[2][1],
                    self.m[2][0] * rhs.m[0][2]
                        + self.m[2][1] * rhs.m[1][2]
                        + self.m[2][2] * rhs.m[2][2],
                ],
            ],            
        }
    }
}

impl<T> Mul<Vector3<T>> for Matrix3<T>
where T: Add<T, Output = T> + Mul<T, Output = T> + Copy {
    type Output = Vector3<T>;

    /// Matrix<T> * Vector3<T>
    /// = Vector3<T>
    /// [r0c0, r0c1, r0c2]   [x]   [x']
    /// [r1c0, r1c1, r1c2] * [y] = [y']
    /// [r2c0, r2c1, r2c2]   [z]   [z']
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [3, 2, 1],
    ///     [1, 3, 2],
    ///     [2, 1, 3]
    /// );
    /// let vec = Vector3::new(1, 3, 2);
    /// assert_eq!(mat * vec, Vector3::new(11, 14, 11));
    /// ```
    fn mul(self, rhs: Vector3<T>) -> Self::Output {
        Vector3(
            rhs.0 * self.m[0][0] + rhs.1 * self.m[0][1] + rhs.2 * self.m[0][2],
            rhs.0 * self.m[1][0] + rhs.1 * self.m[1][1] + rhs.2 * self.m[1][2],
            rhs.0 * self.m[2][0] + rhs.1 * self.m[2][1] + rhs.2 * self.m[2][2],
        )
    }
}

impl<T> Mul<T> for Matrix3<T>
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> * T
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let val = 2;
    /// let expected = Matrix3::new(
    ///     [2, 4, 6],
    ///     [4, 6, 2],
    ///     [6, 2, 4] 
    /// );
    /// assert_eq!(mat * val, expected);
    /// ```
    fn mul(self, other: T) -> Self::Output {
        Matrix3{
            m: [
                [self.m[0][0] * other, self.m[0][1] * other, self.m[0][2] * other],
                [self.m[1][0] * other, self.m[1][1] * other, self.m[1][2] * other],
                [self.m[2][0] * other, self.m[2][1] * other, self.m[2][2] * other],
            ]
        }
    }
}

impl<T> Div for Matrix3<T>
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> / Matrix3<T>
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat0 = Matrix3::new(
    ///     [4, 5, 6],
    ///     [6, 4, 5],
    ///     [5, 6, 4]
    /// );
    /// let mat1 = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let expected = Matrix3::new(
    ///     [1, 0, 0],
    ///     [0, 0, 1],
    ///     [0, 1, 0] 
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

impl<T> Div<T> for Matrix3<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Matrix3<T>;

    /// Matrix3<T> / T
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let val = 2;
    /// let expected = Matrix3::new(
    ///     [0, 1, 1],
    ///     [1, 1, 0],
    ///     [1, 0, 1]
    /// );
    /// assert_eq!(mat / val, expected);
    /// ```
    /// ```should_panic
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// mat / 0;
    /// ```
    fn div(self, rhs: T) -> Self::Output {
        Matrix3 {
            m: [
                [self.m[0][0] / rhs, self.m[0][1] / rhs, self.m[0][2] / rhs],
                [self.m[1][0] / rhs, self.m[1][1] / rhs, self.m[1][2] / rhs],
                [self.m[2][0] / rhs, self.m[2][1] / rhs, self.m[2][2] / rhs],
            ],
        }
    }
}


impl<T> Rem for Matrix3<T>
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    /// Matrix3<T> % Matrix3<T>
    /// = Matrix3<T>
    /// Finds the remainder after an ELEMENT WISE division
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat0 = Matrix3::new(
    ///     [4, 5, 6],
    ///     [6, 4, 5],
    ///     [5, 6, 4]
    /// );
    /// let mat1 = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let expected = Matrix3::new(
    ///     [0, 1, 0],
    ///     [0, 1, 0],
    ///     [2, 0, 0] 
    /// );
    /// assert_eq!(mat0 % mat1, expected);
    /// ```
    fn rem(self, rhs: Self) -> Self::Output {
        Matrix3::new(
            [self.m[0][0] % rhs.m[0][0], self.m[0][1] % rhs.m[0][1], self.m[0][2] % rhs.m[0][2]],
            [self.m[1][0] % rhs.m[1][0], self.m[1][1] % rhs.m[1][1], self.m[1][2] % rhs.m[1][2]],
            [self.m[2][0] % rhs.m[2][0], self.m[2][1] % rhs.m[2][1], self.m[2][2] % rhs.m[2][2]],
        )
    }
}

impl<T> Rem<T> for Matrix3<T> 
where T: Rem<T, Output = T> + Copy {
    type Output = Matrix3<T>;

    /// Matrix3<T> % T
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let val = 2;
    /// let expected = Matrix3::new(
    ///     [1, 0, 1],
    ///     [0, 1, 1],
    ///     [1, 1, 0]
    /// );
    /// assert_eq!(mat % val, expected);
    /// ```
    /// ```should_panic
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// mat % 0;
    /// ```
    fn rem(self, rhs: T) -> Self::Output {
        Matrix3 {
            m: [
                [self.m[0][0] % rhs, self.m[0][1] % rhs, self.m[0][2] % rhs],
                [self.m[1][0] % rhs, self.m[1][1] % rhs, self.m[1][2] % rhs],
                [self.m[2][0] % rhs, self.m[2][1] % rhs, self.m[2][2] % rhs],
            ],
        }
    }
}

impl<T> Neg for Matrix3<T>
where T: Neg<Output = T> + Copy {
    type Output = Self;
    
    /// -Matrix3<T>
    /// = Matrix3<T>
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::matrix3::Matrix3;
    /// let mat = Matrix3::new(
    ///     [1, 2, 3],
    ///     [2, 3, 1],
    ///     [3, 1, 2]
    /// );
    /// let expected = Matrix3::new(
    ///     [-1, -2, -3],
    ///     [-2, -3, -1],
    ///     [-3, -1, -2]
    /// );
    /// assert_eq!(-mat, expected);
    /// ```
    fn neg(self) -> Self::Output {
        Matrix3{
            m: [
                [-self.m[0][0], -self.m[0][1], -self.m[0][2]],
                [-self.m[1][0], -self.m[1][1], -self.m[1][2]],
                [-self.m[2][0], -self.m[2][1], -self.m[2][2]],
            ],
        }
    }
}

impl<T> From<Matrix3<T>> for Matrix4<T>
where T: One + Zero + Copy {
    fn from(mat3: Matrix3<T>) -> Matrix4<T> {
        Matrix4::new(
            [mat3.m[0][0], mat3.m[0][1], mat3.m[0][2], T::ZERO],
            [mat3.m[1][0], mat3.m[1][1], mat3.m[1][2], T::ZERO],
            [mat3.m[2][0], mat3.m[2][1], mat3.m[2][2], T::ZERO],
            [     T::ZERO,      T::ZERO,      T::ZERO,  T::ONE]
        )
    }
}

impl<T> Zero for Matrix3<T> where T: Zero + Copy {
    const ZERO: Self = Matrix3 {
        m: [
            [T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO]
        ]
    };
}

impl<T> One for Matrix3<T> where T: Zero + One {
    /// The identity matrix
    const ONE: Self = Matrix3 {
        m: [
            [T::ONE , T::ZERO, T::ZERO],
            [T::ZERO, T::ONE , T::ZERO],
            [T::ZERO, T::ZERO, T::ONE ]
        ]
    };
}

impl<T> Two for Matrix3<T> where T: Zero + Two {
    /// The identity matrix * 2 
    const TWO: Self = Matrix3 {
        m: [
            [T::TWO , T::ZERO, T::ZERO],
            [T::ZERO, T::TWO , T::ZERO],
            [T::ZERO, T::ZERO, T::TWO ]
        ]
    };
}

impl<T> From<Matrix3<T>> for [[T; 3]; 3]{
    fn from(other: Matrix3<T>) -> [[T; 3]; 3] {
        other.m
    }
}

impl<T> F32Fmt for Matrix3<T> 
where T: F32Fmt + Copy{ 
    type F32Fmt = Matrix3<T::F32Fmt>;
    #[inline]
    fn intoF32Fmt(self) -> Self::F32Fmt {
        let m = &self.m;
        Matrix3{
            m: [
                [m[0][0].intoF32Fmt(), m[0][1].intoF32Fmt(), m[0][2].intoF32Fmt()],
                [m[1][0].intoF32Fmt(), m[1][1].intoF32Fmt(), m[1][2].intoF32Fmt()],
                [m[2][0].intoF32Fmt(), m[2][1].intoF32Fmt(), m[2][2].intoF32Fmt()],
            ]
        }
    }
    #[inline]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { 
        let m = &f32_fmt.m;
        Matrix3{
            m: [
                [T::fromF32Fmt(m[0][0]), T::fromF32Fmt(m[0][1]), T::fromF32Fmt(m[0][2])],
                [T::fromF32Fmt(m[1][0]), T::fromF32Fmt(m[1][1]), T::fromF32Fmt(m[1][2])],
                [T::fromF32Fmt(m[2][0]), T::fromF32Fmt(m[2][1]), T::fromF32Fmt(m[2][2])],
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
        Matrix3{
            m: [
                [self.m[0][0].f32_const_mul(constant), self.m[0][1].f32_const_mul(constant), self.m[0][2].f32_const_mul(constant)],
                [self.m[1][0].f32_const_mul(constant), self.m[1][1].f32_const_mul(constant), self.m[1][2].f32_const_mul(constant)],
                [self.m[2][0].f32_const_mul(constant), self.m[2][1].f32_const_mul(constant), self.m[2][2].f32_const_mul(constant)],
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

impl<T> SignOps for Matrix3<T> {
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