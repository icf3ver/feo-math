//! Constructs for representing rotation in 3D space

use crate::{One, Zero, axes::Axes, linear_algebra::vector3::Vector3};

pub mod quaternion;
pub mod euler;
pub mod rotor;
pub mod axes;

use {
    euler::Euler,
    rotor::Rotor,
    quaternion::Quaternion,

    crate::linear_algebra::matrix3::Matrix3,
};

pub enum RotationTy<T>{
    Quaternion(Quaternion<T>),
    Euler(Euler<T>),
    Rotor(Rotor<T>),
    Axes(Axes<T>)
}

pub trait Rotation<T>:
        Into<Matrix3<T>> +
        From<Matrix3<T>> + 
        Into<Axes<T>> + From<Axes<T>> + 
        Into<Rotor<T>> + From<Rotor<T>> + 
        Into<Quaternion<T>> + From<Quaternion<T>> where T: One + Zero /* does not mean the same for all */{
    const X_AXIS: Vector3<T> = Vector3(T::ONE , T::ZERO, T::ZERO);
    const Y_AXIS: Vector3<T> = Vector3(T::ZERO, T::ONE , T::ZERO);
    const Z_AXIS: Vector3<T> = Vector3(T::ZERO, T::ZERO, T::ONE );

    /// Generate a rotation for an object at a position look at a point.
    /// only uses the X and Y rotation axes to arrive at this result
    /// Z rotation axis
    /// ```text
    ///  ↻
    ///  |  locked Z rotation axis
    ///  | /
    ///  o ---- ⤺
    /// ```
    ///  # Arguments
    /// * `pos` - A `Vector3<f32>` representing the position of the object
    /// * `look_at` - A `Vector3<f32>` representing the position the object is to look at
    fn look_at_xy(pos: Vector3<T>, look_at: Vector3<T>) -> Self;
    
    /// Generate a rotation for an object at a position look at a point.
    /// only uses the X and Z rotation axes to arrive at this result
    /// Z rotation axis
    /// ```text
    ///  locked Y rotation axis
    ///  |  ⤺
    ///  | /
    ///  o ---- ⤺
    /// ```
    ///  # Arguments
    /// * `pos` - A `Vector3<f32>` representing the position of the object
    /// * `look_at` - A `Vector3<f32>` representing the position the object is to look at
    fn look_at_xz(pos: Vector3<T>, look_at: Vector3<T>) -> Self;
    
    /// Generate a rotation for an object at a position look at a point.
    /// only uses the Z and Y rotation axes to arrive at this result
    /// Z rotation axis
    /// ```text
    ///  ↻
    ///  |  ⤺
    ///  | /
    ///  o ---- locked X rotation axis
    /// ```
    ///  # Arguments
    /// * `pos` - A `Vector3<f32>` representing the position of the object
    /// * `look_at` - A `Vector3<f32>` representing the position the object is to look at
    fn look_at_yz(pos: Vector3<T>, look_at: Vector3<T>) -> Self;
    
    /// Generate a rotation for an object at a position look at a point.
    /// only uses the Z and Y rotation axes to arrive at this result
    /// Z rotation axis
    /// ```text
    ///  ↻
    ///  |  ⤺
    ///  | /
    ///  o ---- ⤺
    ///   \
    ///    locked_axis
    /// ```
    ///  # Arguments
    /// * `pos` - A `Vector3<f32>` representing the position of the object
    /// * `look_at` - A `Vector3<f32>` representing the position the object is to look at
    fn look_at_lock(pos: Vector3<T>, look_at: Vector3<T>, locked_axis: Vector3<T>) -> Self;
    
    /// Camera Objects look along the -Z axis so the look_at function for a camera object
    /// needs to be modified a little. flipping the pos and the look_at should do the trick.
    /// TODO: finish comments watch out bc z is reversed and things may get wonky
    /// ```text
    ///    ^ Z
    ///    |
    ///    |
    /// Y• C ------> X
    ///   / \
    ///  /   \
    /// /     \
    /// ```
    /// # Arguments 
    /// * `pos` - A `Vector3<f32>` representing the position of the camera
    /// * `look_at` - A `Vector3<f32>` representing the position the camera is to look at
    fn camera_look_at_xy(pos: Vector3<T>, look_at: Vector3<T>) -> Self;
    fn camera_look_at_xz(pos: Vector3<T>, look_at: Vector3<T>) -> Self;
    fn camera_look_at_yz(pos: Vector3<T>, look_at: Vector3<T>) -> Self;
    fn camera_look_at_lock(pos: Vector3<T>, look_at: Vector3<T>, locked_axis: Vector3<T>) -> Self;
}


