//! You can represent a rotation through the rotation of axes.
//! 
//! TODO: explain
//! 
use std::{fmt, ops::{Add, Div, Mul, Neg, Sub}};
use crate::{F32Fmt, One, Two, Zero, linear_algebra::{matrix3::Matrix3, vector3::Vector3}, rotation::quaternion::Quaternion};

use super::{euler::Euler, rotor::Rotor};

pub static NORMAL_X_AXIS: Vector3<f32>  = Vector3( 1.0,  0.0,  0.0);
pub static NORMAL_Y_AXIS: Vector3<f32>  = Vector3( 0.0,  1.0,  0.0);
pub static NORMAL_Z_AXIS: Vector3<f32>  = Vector3( 0.0,  0.0,  1.0);

pub static VULKANO_Y_AXIS: Vector3<f32> = Vector3( 0.0, -1.0,  0.0);

#[derive(Clone, Copy)]
pub struct Axes<T>{
    pub x: Vector3<T>, // left
    pub y: Vector3<T>, // up
    pub z: Vector3<T>  // frd or forward
}

impl<T> Axes<T> {
    /// Returns a set of axes calculated using the passed in x and z axes
    /// # Arguments
    /// * `left` - The x axis represented by a Vector3
    /// * `frd` - The z axis represented by a Vector3
    pub fn from_left_frd(left: Vector3<T>, frd: Vector3<T>) -> Self 
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Copy{
        let up = Vector3::cross_product(frd, left);

        let rot_mat =  Matrix3{
            m: [
                left.into(),
                up.into(),
                frd.into()
            ]
        }.transpose();
        
        rot_mat.into()
    }

    /// Returns a set of axes calculated using the passed in y and z axes
    /// # Arguments
    /// * `up` - The y axis represented by a Vector3
    /// * `frd` - The z axis represented by a Vector3
    pub fn from_up_frd(up: Vector3<T>, frd: Vector3<T>) -> Self 
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Copy {
        let left = Vector3::cross_product(up, frd);

        let rot_mat = Matrix3{
            m: [
                left.into(),
                up.into(),
                frd.into()
            ]
        }.transpose(); // not needed TODO
        
        rot_mat.into()
    }

    /// Returns a set of axes calculated using the passed in x and y axes
    /// # Arguments
    /// * `left` - The x axis represented by a Vector3
    /// * `up` - The y axis represented by a Vector3
    pub fn from_left_up(left: Vector3<T>, up: Vector3<T>) -> Self 
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Copy {
        let frd = Vector3::cross_product(left, up);

        let rot_mat = Matrix3{
            m: [
                left.into(),
                up.into(),
                frd.into()
            ]
        }.transpose();
        
        rot_mat.into()
    }
}

impl<T> fmt::Debug for Axes<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Axes: 
                 \n   x: {:?} 
                 \n   y: {:?} 
                 \n   z: {:?}", self.x, self.y, self.z)
    }
}

impl<T> From<Quaternion<T>> for Axes<T>
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + F32Fmt + fmt::Debug + One + Two + Zero + Copy{
    fn from(q: Quaternion<T>) -> Self {
        Axes::from(Into::<Matrix3<T>>::into(q))
    }
}

impl<T> From<Euler<T>> for Axes<T> 
where T: Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Add<T, Output = T> + Copy {
    fn from(e: Euler<T>) -> Self {
        Axes::from(Into::<Matrix3<T>>::into(e))
    }
}

impl<T> From<Rotor<T>> for Axes<T> 
where T: Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Add<T, Output = T> + Copy {
    fn from(r: Rotor<T>) -> Self {
        Axes::from(Into::<Matrix3<T>>::into(r))
    }
}

impl<T> From<Matrix3<T>> for Axes<T>
where T: Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Add<T, Output = T> + Copy {
    
        // so the transposed rotation matrix is:
        // [  left ]
        // [  up   ] 
        // [  frd  ]

    fn from(mat: Matrix3<T>) -> Self {
        let mat = mat.transpose();
        Axes {
            x: Vector3::from(mat.m[0]).normalize(None),
            y: Vector3::from(mat.m[1]).normalize(None),
            z: Vector3::from(mat.m[2]).normalize(None),
        }
    }
}

impl<T> From<Axes<T>> for Matrix3<T> 
where T: Copy {
    fn from(other: Axes<T>) -> Matrix3<T> {
        Matrix3::new(
            other.x.into(),
            other.y.into(),
            other.z.into(),
        ).transpose()
    }
}