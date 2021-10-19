use core::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
/// call mul_var

#[derive(Copy, Clone)]
pub struct WithVar<T, C = char>{
    pub val: T,
    pub var: C,
}

impl<T> fmt::Debug for WithVar<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}{}", self.val, self.var)
    }
}

impl<T, C> Add for WithVar<T, C>
where T: Add<T, Output = T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        WithVar{
            val: self.val + rhs.val,
            var: self.var,
        }
    }
}

impl<T, C> Sub for WithVar<T, C>
where T: Sub<T, Output = T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        WithVar{
            val: self.val - rhs.val,
            var: self.var,
        }
    }
}
impl<T, C> Mul for WithVar<T, C>
where T: Mul<T, Output = T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        WithVar{
            val: self.val * rhs.val,
            var: self.var,
        }
    }
}

impl<T, C> Div for WithVar<T, C>
where T: Div<T, Output = T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        WithVar{
            val: self.val / rhs.val,
            var: self.var,
        }
    }
}
impl<T, C> Neg for WithVar<T, C>
where T: Neg<Output = T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        WithVar{
            val: -self.val,
            var: self.var,
        }
    }
}

/// ...
/// ...

impl<T, C> Mul<T> for WithVar<T, C>
where T: Mul<T, Output = T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        WithVar{
            val: self.val * rhs,
            var: self.var,
        }
    }
}

impl<T, C> Div<T> for WithVar<T, C>
where T: Div<T, Output = T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        WithVar{
            val: self.val / rhs,
            var: self.var,
        }
    }
}