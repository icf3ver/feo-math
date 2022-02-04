//! Construct that represents euler rotations.
//! 
//! Very simple but flawed.
//! TODO
//! 
use std::{fmt, ops::{Add, Div, Mul, Neg, Rem, Sub}};

use crate::{Construct, F32Fmt, One, SignOps, /*Sum,*/ Two, Zero, axes::Axes, linear_algebra::{matrix3::Matrix3, vector3::Vector3}};

use super::{Rotation, quaternion::Quaternion, rotor::Rotor};

/// Rotations around the x y and z axes in radians radian = value/PI do not include PI
#[derive(Copy, Clone)]
pub struct Euler<T>(pub T, pub T, pub T);

impl<T> PartialEq for Euler<T>
where Self: Rem<T, Output = Self>, T: Add<T, Output = T> + One + PartialEq + Copy {
    fn eq(&self, other: &Self) -> bool {
        let two = T::ONE + T::ONE;
        let spa0 = *self % two;
        let spa1 = *other % two;
        spa0.0 == spa1.0 && spa0.1 == spa1.1 && spa0.2 == spa1.2
    }
}

impl<T> Euler<T>{
    pub fn new(x_rot: T, y_rot: T, z_rot: T) -> Self {
        Euler(x_rot, y_rot, z_rot)
    }
}

impl<T> Construct<T> for Euler<T> where T: Construct<T> {}
impl<T> Rotation<T> for Euler<T> where T: Construct<T> {
    fn look_at_xy(_pos: Vector3<T>, _look_at: Vector3<T>) -> Self {
        todo!()
    }

    fn look_at_xz(_pos: Vector3<T>, _look_at: Vector3<T>) -> Self {
        todo!()
    }

    fn look_at_yz(_pos: Vector3<T>, _look_at: Vector3<T>) -> Self {
        todo!()
    }

    fn look_at_lock(_pos: Vector3<T>, _look_at: Vector3<T>, _locked_axis: Vector3<T>) -> Self {
        todo!()
    }

    fn camera_look_at_xy(_pos: Vector3<T>, _look_at: Vector3<T>) -> Self {
        todo!()
    }

    fn camera_look_at_xz(_pos: Vector3<T>, _look_at: Vector3<T>) -> Self {
        todo!()
    }

    fn camera_look_at_yz(_pos: Vector3<T>, _look_at: Vector3<T>) -> Self {
        todo!()
    }

    fn camera_look_at_lock(_pos: Vector3<T>, _look_at: Vector3<T>, _locked_axis: Vector3<T>) -> Self {
        todo!()
    }
}

impl<T> fmt::Debug for Euler<T> 
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Euler({:?}, {:?}, {:?})", self.0, self.1, self.2)
    }
}

impl<T> From<Euler<T>> for Matrix3<T> {
    fn from(_: Euler<T>) -> Matrix3<T> {
        todo!()
    }
}

impl<T> From<Matrix3<T>> for Euler<T> {
    fn from(_: Matrix3<T>) -> Self {
        todo!()
    }
}

impl<T> From<Axes<T>> for Euler<T> 
where T: Copy {
    fn from(a: Axes<T>) -> Self {
        Into::<Matrix3<T>>::into(a).into()
    }
}

impl<T> From<Rotor<T>> for Euler<T> 
where T: Construct<T> + Copy {
    fn from(r: Rotor<T>) -> Self {
        Into::<Matrix3<T>>::into(r).into()
    }
}

impl<T> From<Quaternion<T>> for Euler<T> 
where T: Construct<T> + Copy {
    fn from(q: Quaternion<T>) -> Self {
        Into::<Matrix3<T>>::into(q).into()
    }
}

impl<T> Add for Euler<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: Euler<T>) -> Self::Output {
        Euler(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2)
    }
}

impl<T> Add<T> for Euler<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        Euler(self.0 + rhs, self.1 + rhs, self.2 + rhs)
    }
}

impl<T> Sub for Euler<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: Euler<T>) -> Self::Output {
        Euler(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2)
    }
}

impl<T> Sub<T> for Euler<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        Euler(self.0 - rhs, self.1 - rhs, self.2 - rhs)
    }
}

impl<T> Mul for Euler<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: Euler<T>) -> Self::Output {
        Euler(self.0 * rhs.0, self.1 * rhs.1, self.2 * rhs.2)
    }
}

impl<T> Mul<T> for Euler<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Euler(self.0 * rhs, self.1 * rhs, self.2 * rhs)
    }
}

impl<T> Div for Euler<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: Euler<T>) -> Self::Output {
        Euler(self.0 / rhs.0, self.1 / rhs.1, self.2 / rhs.2)
    }
}

impl<T> Div<T> for Euler<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Euler(self.0 / rhs, self.1 / rhs, self.2 / rhs)
    }
}

impl<T> Rem for Euler<T> 
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    fn rem(self, rhs: Euler<T>) -> Self::Output {
        Euler(self.0 % rhs.0, self.1 % rhs.1, self.2 % rhs.2)
    }
}

impl<T> Rem<T> for Euler<T> 
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    fn rem(self, rhs: T) -> Self::Output {
        Euler(self.0 % rhs, self.1 % rhs, self.2 % rhs)
    }
}

impl<T> Neg for Euler<T> 
where T: Neg<Output = T> + Copy {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Euler(-self.0, -self.1, -self.2)
    }
}

impl<T> Zero for Euler<T> where T: Zero { const ZERO: Self = Euler(T::ZERO, T::ZERO, T::ZERO); }
impl<T> One for Euler<T> where T: One { const ONE: Self = Euler(T::ONE, T::ONE, T::ONE); }
impl<T> Two for Euler<T> where T: Two { const TWO: Self = Euler(T::TWO, T::TWO, T::TWO); }

impl<T: F32Fmt> F32Fmt for Euler<T> { 
    type F32Fmt = Euler<T::F32Fmt>;
    #[inline]
    fn intoF32Fmt(self) -> Self::F32Fmt {
        Euler(self.0.intoF32Fmt(), self.1.intoF32Fmt(), self.2.intoF32Fmt())
    }
    #[inline]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self {
        Euler(T::fromF32Fmt(f32_fmt.0), T::fromF32Fmt(f32_fmt.1), T::fromF32Fmt(f32_fmt.2))
    }

    fn sqrt(self) -> Self {
        Euler(self.0.sqrt(), self.1.sqrt(), self.2.sqrt())
    }

    fn cbrt(self) -> Self {
        Euler(self.0.cbrt(), self.1.cbrt(), self.2.cbrt())
    }

    fn f32_const_mul(self, constant: f32) -> Self {
        Euler(self.0.f32_const_mul(constant), self.1.f32_const_mul(constant), self.2.f32_const_mul(constant))
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

impl<T> SignOps for Euler<T> {
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