use std::ops::Rem;

use crate::{F32Fmt, One, SignOps, /*Sum,*/ Two, Zero};

use {
    super::Vector,
    crate::Construct,
    std::{
        fmt, 
        ops::{Add, Div, Sub, Mul, Neg}
    }
};

/// A 3D vector
#[derive(PartialEq, Clone, Copy)]
pub struct Vector3<T>(pub T, pub T, pub T);

impl<T> fmt::Debug for Vector3<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:?}, {:?}, {:?})", self.0, self.1, self.2)
    }
}

impl<T> Vector3<T> {
    /// A clean way of making a Vector3 <x, y, z>
    /// # Arguments
    /// * `x` - The x direction scalar .0
    /// * `y` - The y direction scalar .1
    /// * `z` - The z direction scalar .2
    pub fn new(x: T, y: T, z: T) -> Self
    where T: Construct<T> {
        Vector3(x, y, z)
    }

    /// Returns the Vector3<T> normalized to a given length
    /// # Arguments
    /// * `length` - The length to normalize to. Default is 1
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// let vec = Vector3::new(2, 0, 0);
    /// assert_eq!(vec.normalize(Some(5)), Vector3::new(5, 0, 0));
    /// assert_eq!(vec.normalize(None), Vector3::new(1, 0, 0));
    /// ```
    pub fn normalize(&self, length: Option<T>) -> Vector3<T>
    where T: Mul<T, Output = T> + Div<T, Output = T> + Add<T, Output = T> + F32Fmt + Copy {
        let s = self.norm();
        match length {
            Some(len) => Vector3(self.0 / s * len, self.1 / s * len, self.2 / s * len),
            None => Vector3(self.0 / s, self.1 / s, self.2 / s),
        }
    }

    /// Returns the magnitude, norm, or length of the vector denoted  `||v||`
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// let vec = Vector3::new(2, 0, 0);
    /// assert_eq!(vec.norm(), 2);
    /// ```
    #[inline]
    pub fn norm(&self) -> T
    where T: Add<T, Output = T> + Mul<T, Output = T> + F32Fmt + Copy {
        (self.0 * self.0 + self.1 * self.1 + self.2 * self.2).sqrt()
    }

    /// Returns the unit vector of the vector denoted `Uv`
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// let v = Vector3::new(1.0, 2.0, 3.0);
    /// let Uv = v.unit_vector();
    /// assert_eq!(Uv.norm(), 1.0);
    /// ```
    #[inline]
    pub fn unit_vector(&self) -> Self
    where Self: Div<T, Output = Self>, T: Add<T, Output = T> + Mul<T, Output = T> + F32Fmt + Copy {
        *self / self.norm()
    }

    /// Returns the cross product of 2 Vector4s denoted `A x B`
    /// A trick to know which direction the cross product vector
    /// will be pointing is the right hand rule
    /// ```text
    /// side view
    ///          ^\ a x b (pointing up)
    /// a _______|| _________ arm (you)
    ///       ///_  
    ///      | \__  
    ///      |
    ///      b (pointing left)
    /// ```
    /// where your index finger is vector a pointing forwards along z
    /// your middle finger is vector b pointing left along x
    /// and your thumb is the cross product a x b pointing up
    /// If you want to think in standard camera space terms use 
    /// the left hand rule and note that your index finger is pointing 
    /// along -z
    /// # Arguments
    /// * `a` - The first Vector4
    /// * `b` - The second Vector4
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// let a = Vector3(0, 0, -1);
    /// let b = Vector3(-1, 0, 0);
    /// assert_eq!(Vector3::cross_product(a, b), Vector3(0, 1, 0));
    /// ```
    pub fn cross_product(a: Self, b: Self) -> Self
    where T: Mul<T, Output = T> + Sub<T, Output = T> + Copy {
        Vector3(
            (a.1 * b.2) - (a.2 * b.1),
            (a.2 * b.0) - (a.0 * b.2),
            (a.0 * b.1) - (a.1 * b.0)
        )
    }

    /// Returns the dot product of 2 Vector4s denoted `A . B`
    /// # Arguments
    /// * `a` - The first Vector4
    /// * `b` - The second Vector4
    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// let a = Vector3(0, 0, -1);
    /// let b = Vector3(-1, 0, 0);
    /// assert_eq!(Vector3::dot_product(a, b), 0);
    /// ```
    pub fn dot_product(a: Self, b: Self) -> T
    where T: Mul<T, Output = T> + Add<T, Output = T> + Copy {
        (a.0 * b.0) + (a.1 * b.1) + (a.2 * b.2)
    }
}

impl<T> Construct<T> for Vector3<T> where T: Construct<T> {}
impl<T> Vector<T> for Vector3<T> where T: Construct<T> {}

impl<T> Add for Vector3<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    /// Vector
    /// # Examples
    /// ```
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// let vec = Vector3::new(1, 2, 3);
    /// let expected = Vector3::new(2, 4, 6);
    /// assert_eq!(vec + vec, expected);
    /// ```
    fn add (self, rhs: Self) -> Self::Output{
        Vector3(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2)
    }
}

impl<T> Add<T> for Vector3<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    /// Element wise addition
    /// # Examples
    /// ```
    /// use feo_math::linear_algebra::vector3::Vector3;
    /// assert_eq!(Vector3(1, 3, 2) + 2, Vector3(3, 5, 4));
    /// ```
    fn add (self, rhs: T) -> Self::Output{
        Vector3(self.0 + rhs, self.1 + rhs, self.2 + rhs)
    }
}

impl<T> Sub for Vector3<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub (self, rhs: Self) -> Self::Output{
        Vector3(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2)
    }
}

impl<T> Sub<T> for Vector3<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub (self, rhs: T) -> Self::Output{
        Vector3(self.0 - rhs, self.1 - rhs, self.2 - rhs)
    }
}

impl<T> Mul for Vector3<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul (self, rhs: Self) -> Self::Output{
        Vector3(self.0 * rhs.0, self.1 * rhs.1, self.2 * rhs.2)
    }
}

impl<T> Mul<T> for Vector3<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul (self, rhs: T) -> Self::Output{
        Vector3(self.0 * rhs, self.1 * rhs, self.2 * rhs)
    }
}

impl<T> Div for Vector3<T>
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div (self, rhs: Self) -> Self::Output{
        Vector3(self.0 / rhs.0, self.1 / rhs.1, self.2 / rhs.2)
    }
}

impl<T> Div<T> for Vector3<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div (self, rhs: T) -> Self::Output {
        Vector3(self.0 / rhs, self.1 / rhs, self.2 / rhs)
    }
}

impl<T> Rem for Vector3<T>
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    fn rem (self, rhs: Self) -> Self::Output {
        Vector3(self.0 % rhs.0, self.1 % rhs.1, self.2 % rhs.2)
    }
}

impl<T> Rem<T> for Vector3<T> 
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    /// # Examples
    /// ```rust
    /// use feo_math::linear_algebra::vector3;
    /// ```
    fn rem (self, rhs: T) -> Self::Output {
        Vector3(self.0 % rhs, self.1 % rhs, self.2 % rhs)
    }
}

impl<T> Neg for Vector3<T> 
where T: Neg<Output = T> + Copy {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Vector3(-self.0, -self.1, -self.2)
    }
}

impl<T> Zero for Vector3<T> where T: Zero { const ZERO: Self = Vector3(T::ZERO, T::ZERO, T::ZERO); }
impl<T> One for Vector3<T> where T: One { const ONE: Self = Vector3(T::ONE, T::ONE, T::ONE); }
impl<T> Two for Vector3<T> where T: Two { const TWO: Self = Vector3(T::TWO, T::TWO, T::TWO); }

impl<T> From<Vector3<T>> for [T; 3] 
where T: Copy {
    fn from(other: Vector3<T>) -> [T; 3]{
        [other.0, other.1, other.2]
    }
}

impl<T> From<[T; 3]> for Vector3<T> 
where T: Copy{
    fn from(other: [T; 3]) -> Self {
        Vector3(other[0], other[1], other[2])
    }
}

impl<T> F32Fmt for Vector3<T> 
where   T: F32Fmt + Copy { 
    type F32Fmt = Vector3<T::F32Fmt>;
    #[inline]
    fn intoF32Fmt(self) -> Self::F32Fmt {
        Vector3(self.0.intoF32Fmt(), self.1.intoF32Fmt(), self.2.intoF32Fmt())
    }
    #[inline]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self {
        let vec = &f32_fmt;
        Vector3(T::fromF32Fmt(vec.0), T::fromF32Fmt(vec.1), T::fromF32Fmt(f32_fmt.2))
    }

    /// Element wise square root
    fn sqrt(self) -> Self {
       Vector3(self.0.sqrt(), self.1.sqrt(), self.2.sqrt())
    }

    fn cbrt(self) -> Self { 
        Vector3(self.0.cbrt(), self.1.cbrt(), self.2.cbrt())
    }

    fn f32_const_mul(self, constant: f32) -> Self {
        Vector3(self.0.f32_const_mul(constant), self.1.f32_const_mul(constant), self.2.f32_const_mul(constant))
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
impl<T> SignOps for Vector3<T> {
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