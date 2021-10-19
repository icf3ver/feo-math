use std::{fmt, ops::{Add, Div, Mul, Neg, Rem, Sub}};

use crate::{Construct, F32Fmt, One, SignOps, /*Sum,*/ Two, Zero};

use super::{Vector, vector3::Vector3};

/// A 4D vector
#[derive(PartialEq, Clone, Copy)]
pub struct Vector4<T>(pub T, pub T, pub T, pub T);

impl<T> Vector4<T> where T: Copy{
    /// Returns a new Vector4
    /// # Arguments
    /// * `x` - The x direction scalar .0
    /// * `y` - The y direction scalar .1
    /// * `z` - The z direction scalar .2
    /// * `w` - The w direction scalar .3
    pub fn new(x: T, y: T, z: T, w: T) -> Self {
        Vector4(x, y, z, w)
    }

    /// Returns a new Vector4 created by adding a w value to a Vector3
    /// # Arguments
    /// * `vec3` - The 3D Vector that the w direction scaler is added onto
    /// * `w` - The w direction scalar .3
    pub fn from_vector3(vec3: Vector3<T>, w: T) -> Self {
        Vector4(vec3.0, vec3.1, vec3.2, w)
    }

    pub fn empty(zero: T) -> Self {
        Vector4(zero, zero, zero, zero)
    }

    pub fn normalize(&self, length: Option<T>) -> Vector4<T>
    where T: Mul<T, Output = T> + Div<T, Output = T> + Add<T, Output = T> + F32Fmt {
        let size: T = self.norm();
        match length{
            Some(len) => Vector4((self.0 / size) * len, (self.1 / size) * len, (self.2 / size) * len, (self.3 / size) * len),
            None => Vector4(self.0 / size, self.1 / size, self.2 / size, self.3 / size),
        }
    }

    /// Returns the magnitude or length of the Vector4 denoted `||q||`
    ///
    pub fn norm(&self) -> T
    where T: Add<T, Output = T> + Mul<T, Output = T> + F32Fmt + Copy {
        (self.0 * self.0 + self.1 * self.1 + self.2 * self.2 + self.3 * self.3).sqrt()
    }
}

impl<T> Construct<T> for Vector4<T> where T: Construct<T> {}
impl<T> Vector<T> for Vector4<T> where T: Construct<T> {/* none */}

impl<T> fmt::Debug for Vector4<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:?}, {:?}, {:?}, {:?})", self.0, self.1, self.2, self.3)
    }
}

impl<T> Add for Vector4<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add (self, rhs: Self) -> Self::Output {
        Vector4(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2, self.3 + rhs.3)
    }
}

impl<T> Add<T> for Vector4<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add (self, rhs: T) -> Self::Output {
        Vector4(self.0 + rhs, self.1 + rhs, self.2 + rhs, self.3 + rhs)
    }
}

impl<T> Sub for Vector4<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub (self, rhs: Self) -> Self::Output {
        Vector4(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2, self.3 - rhs.3)
    }
}

impl<T> Sub<T> for Vector4<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub (self, rhs: T) -> Self::Output {
        Vector4(self.0 - rhs, self.1 - rhs, self.2 - rhs, self.3 - rhs)
    }
}

impl<T> Mul for Vector4<T> where T: Mul<T, Output = T> + Copy{
    type Output = Self;

    fn mul (self, rhs: Self) -> Self::Output{
        Vector4(self.0 * rhs.0, self.1 * rhs.1, self.2 * rhs.2, self.3 * rhs.3)
    }
}

impl<T> Mul<T> for Vector4<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul (self, rhs: T) -> Self::Output{
        Vector4(self.0 * rhs, self.1 * rhs, self.2 * rhs, self.3 * rhs)
    }
}

impl<T> Div for Vector4<T>
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div (self, rhs: Self) -> Self::Output{
        Vector4(self.0 / rhs.0, self.1 / rhs.1, self.2 / rhs.2, self.3 / rhs.3)
    }
}

impl<T> Div<T> for Vector4<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div (self, rhs: T) -> Self::Output{
        Vector4(self.0 / rhs, self.1 / rhs, self.2 / rhs, self.3 / rhs)
    }
}

impl<T> Rem for Vector4<T>
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    /// Vector4<T> % Vector4<T>
    /// = Vector4<T>
    /// # Examples
    /// ```
    /// use feo_math::linear_algebra::vector4::Vector4;
    /// 
    /// ```
    fn rem (self, rhs: Self) -> Self::Output{
        Vector4(self.0 % rhs.0, self.1 % rhs.1, self.2 % rhs.2, self.3 % rhs.3)
    }
}

impl<T> Rem<T> for Vector4<T> 
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    fn rem (self, rhs: T) -> Self::Output{
        Vector4(self.0 % rhs, self.1 % rhs, self.2 % rhs, self.3 % rhs)
    }
}

impl<T> Neg for Vector4<T> 
where T: Neg<Output = T> + Copy {
    type Output = Self;

    fn neg(self) -> Self::Output{
        Vector4(-self.0, -self.1, -self.2, -self.3)
    }
}

impl<T> Zero for Vector4<T> 
where T: Zero {
    const ZERO: Self = Vector4(T::ZERO, T::ZERO, T::ZERO, T::ZERO);
}
impl<T> One for Vector4<T> 
where T: One {
    const ONE: Self = Vector4(T::ONE, T::ONE, T::ONE, T::ONE);
}
impl<T> Two for Vector4<T> 
where T: Two {
    const TWO: Self = Vector4(T::TWO, T::TWO, T::TWO, T::TWO);
}

impl<T> From<Vector4<T>> for [T; 4]{
    fn from(other: Vector4<T>) -> [T; 4]{
        [other.0, other.1, other.2, other.3]
    }
}

impl<T> From<Vector4<T>> for Vector3<T> where T: Copy{
    fn from(other: Vector4<T>) -> Vector3<T>{
        Vector3(other.0, other.1, other.2)
    }
}

impl<T> From<[T; 4]> for Vector4<T> 
where T: Copy{
    fn from(other: [T; 4]) -> Self {
        Vector4(other[0], other[1], other[2], other[3])
    }
}

impl<T: F32Fmt> F32Fmt for Vector4<T> 
where T: Mul<T, Output = T> { 
    type F32Fmt = Vector4<T::F32Fmt>;
    #[inline]
    fn intoF32Fmt(self) -> Self::F32Fmt {
        Vector4(self.0.intoF32Fmt(), self.1.intoF32Fmt(), self.2.intoF32Fmt(), self.3.intoF32Fmt())
    } 
    #[inline]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self {
        let vec = &f32_fmt;
        Vector4(T::fromF32Fmt(vec.0), T::fromF32Fmt(vec.1), T::fromF32Fmt(f32_fmt.2), T::fromF32Fmt(f32_fmt.3))
    }

    fn sqrt(self) -> Self {
        Vector4(self.0.sqrt(), self.1.sqrt(), self.2.sqrt(), self.3.sqrt())
    }

    fn cbrt(self) -> Self {
        Vector4(self.0.cbrt(), self.1.cbrt(), self.2.cbrt(), self.3.cbrt())
    }

    fn f32_const_mul(self, constant: f32) -> Self {
        Vector4(self.0.f32_const_mul(constant), self.1.f32_const_mul(constant), self.2.f32_const_mul(constant), self.3.f32_const_mul(constant))
    }

    fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        Vector4(self.0.sin_mul(mul_by.0), self.1.sin_mul(mul_by.1), self.2.sin_mul(mul_by.2), self.3.sin_mul(mul_by.3))
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
impl<T> SignOps for Vector4<T> {
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