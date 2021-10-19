use crate::{Construct, One, Two, Zero, axes::Axes, linear_algebra::{matrix3::Matrix3, vector3::Vector3}};

use super::{Rotation, euler::Euler, quaternion::Quaternion};

/// Uses the same principles as quaternions to arrive at the same result in a fashion that is easier to understand
pub struct Rotor<T>(pub T, pub T, pub T, pub T);

impl<T> Rotor<T> {
    pub fn new(x_rot: T, y_rot: T, z_rot: T, w_rot: T) -> Self {
        Rotor(x_rot, y_rot, z_rot, w_rot)
    }
}

impl<T> Rotation<T> for Rotor<T> where T: Construct<T> {
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

    fn camera_look_at_xy(pos: Vector3<T>, look_at: Vector3<T>) -> Self { Rotor::look_at_xy(look_at, pos) }

    fn camera_look_at_xz(pos: Vector3<T>, look_at: Vector3<T>) -> Self { Rotor::look_at_xz(look_at, pos) }

    fn camera_look_at_yz(pos: Vector3<T>, look_at: Vector3<T>) -> Self { Rotor::look_at_yz(look_at, pos) }

    fn camera_look_at_lock(pos: Vector3<T>, look_at: Vector3<T>, locked_axis: Vector3<T>) -> Self {
        Rotor::look_at_lock(look_at, pos, locked_axis)
    }
}

impl<T> From<Rotor<T>> for Matrix3<T> {
    fn from(_: Rotor<T>) -> Matrix3<T> {
        todo!()
    }
}

impl<T> From<Matrix3<T>> for Rotor<T> {
    fn from(_: Matrix3<T>) -> Self {
        todo!()
    }
}

impl<T> From<Axes<T>> for Rotor<T>
where T: Copy {
    fn from(other: Axes<T>) -> Self {
        Into::<Matrix3<T>>::into(other).into()
    }
}

impl<T> From<Euler<T>> for Rotor<T> 
where T: Construct<T> + Copy {
    fn from(e: Euler<T>) -> Self {
        Into::<Matrix3<T>>::into(e).into()
    }
}

impl<T> From<Quaternion<T>> for Rotor<T> 
where T: Construct<T> + Copy {
    fn from(q: Quaternion<T>) -> Self {
        Into::<Matrix3<T>>::into(q).into()
    }
}
impl<T> Zero for Rotor<T> where T: Zero { const ZERO: Self = Rotor(T::ZERO, T::ZERO, T::ZERO, T::ZERO); }
impl<T> One for Rotor<T> where T: One + Zero { const ONE: Self = Rotor(T::ZERO, T::ZERO, T::ZERO, T::ONE); }
impl<T> Two for Rotor<T> where T: Two + Zero { const TWO: Self = Rotor(T::ZERO, T::ZERO, T::ZERO, T::TWO); }