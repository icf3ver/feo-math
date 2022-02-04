//! Space a construct for storing everything about how an object is to be transformed.
//! 
//! A space is an area in which objects, such as spaces or the vertices of a mesh, can exist. 
//! An issue with this visualization is however that the area itself is unlimited in size. If 
//! a space is translated, rotated, or scaled the objects within it do not change their 
//! position, rotation, or scaling relative to the space. Their position relative to the space
//!  containing that space however will change. 
//!
//! > Spaces can also be joined to represent conversions between spaces.
//!
//! Standard examples of spaces include:
//!  + Worldspace : The space in which all objects exist without their own space and any spaces containing their own space.
//!  + Viewspace : The space the viewer perceives viewed objects in
//!  + Cameraspace : The position relative to the camera. NOT the viewspace.
//!  + Objectspace : every game-object has its own space
//!
//! 
use std::ops::{Add, Sub};

use crate::{Zero, linear_algebra::{
        matrix3::Matrix3, 
        matrix4::Matrix4, 
        vector3::Vector3
    }, rotation::{Rotation, quaternion::Quaternion}};

#[derive(Debug, Clone, Copy)]
pub struct Space {
    pub center: Vector3<f32>,
    pub rotation: Quaternion<f32>,
    pub scale_factor: Matrix3<f32>
}
// todo move to src no utils should exist
impl Space {
    #[allow(clippy::redundant_closure)] // Linter error: .unwrap_or_else() wants FnOnce<_> closure
    pub fn new(
            center: Option<Vector3<f32>>, 
            rotation: Option<Quaternion<f32>>, 
            scale_factor: Option<Vector3<f32>>) -> Space {
        let vec = scale_factor.unwrap_or_else(|| Vector3::new(1.0, 1.0, 1.0));
        Space{
            rotation: rotation.unwrap_or_else(|| Quaternion::f32_identity()),
            center: center.unwrap_or(Vector3::ZERO),
            scale_factor: Matrix3::new(
                [vec.0, 0.0, 0.0],
                [0.0, vec.1, 0.0],
                [0.0, 0.0, vec.2],
            )
        }
    }

    pub fn identity() -> Space {
        Space{
            rotation: Quaternion::f32_identity(),
            center: Vector3::ZERO,
            scale_factor: Matrix3::identity()
        }
    }

    // get a space representing both self and other
    // self is in other
    // move self to other
    // the resulting matrix can represent the points in other's superspace (kind of)
    // you see you now have the representation in others superspace but without the rotation
    // that is only applied when you then go to the superspace and join it to
    
    // Think of a space like a 2D plane that can be moved anywhere and rotated around its center and scaled
    // The joining of two spaces basically means that now the space is moved to where it should be in the 
    // superspace of the other matrix it is being joined with and it has the rotation of other bc:
    // __ no rotation
    // \_ center rotated
    // \
    //  \ points in space rotated to match change
    // Because the other space has some position within its superspace that is also added to the position
    pub fn join(&self, other: Space) -> Space {
        // Diagonalize scl mat
        let sf = self.scale_factor.transpose().m;

        let left = other.rotation * Quaternion::new_vector_real(Vector3::from(sf[0]), 0.0) * other.rotation.reciprocal();
        let up = other.rotation * Quaternion::new_vector_real(Vector3::from(sf[1]), 0.0) * other.rotation.reciprocal();
        let frd = other.rotation * Quaternion::new_vector_real(Vector3::from(sf[2]), 0.0) * other.rotation.reciprocal();

        let scl_mat = Matrix3::new(
            [left.0, up.0, frd.0],
            [left.1, up.1, frd.1],
            [left.2, up.2, frd.2]
        );

        let rotated_center = other.rotation * Quaternion::new_vector_real(self.center, 0.0) * other.rotation.reciprocal();

        Space{
            rotation: other.rotation * self.rotation,
            center: other.center + Vector3(rotated_center.0, rotated_center.1, rotated_center.2),
            scale_factor: Matrix3::from(other.rotation).transpose() * (scl_mat * other.scale_factor)
        }
    }
    
    // self in other
    // move other to self
    // Think of a space like a 2D plane that can be moved anywhere and rotated around its center and scaled
    // The joining in reverse two spaces basically means that now the superspace is moved space is moved to where it should be in the 
    // superspace of the other matrix it is being joined with and it has the rotation of other bc:
    // __ no rotation
    // \_ center rotated
    // \
    //  \ points in space rotated to match change
    // Because the other space has some position within its superspace that is also added to the position
    pub fn join_reverse(&self, other: Space) -> Space {
        let sf = other.scale_factor.transpose().m;

        let left = self.rotation.reciprocal() * Quaternion::new_vector_real(Vector3::from(sf[0]), 0.0) * self.rotation;
        let up = self.rotation.reciprocal() * Quaternion::new_vector_real(Vector3::from(sf[1]), 0.0) * self.rotation;
        let frd = self.rotation.reciprocal() * Quaternion::new_vector_real(Vector3::from(sf[2]), 0.0) * self.rotation;

        let scl_mat = Matrix3::new(
            [left.0, up.0, frd.0],
            [left.1, up.1, frd.1],
            [left.2, up.2, frd.2]
        );

        let rotated_center = self.rotation.reciprocal() * Quaternion::new_vector_real(other.center - self.center, 0.0) * self.rotation;
        Space{
            rotation: self.rotation.reciprocal() * other.rotation,
            center: Vector3(rotated_center.0, rotated_center.1, rotated_center.2),
            scale_factor: Matrix3::from(self.rotation) * (scl_mat * self.scale_factor.inverse())
        }
    }

    pub fn build(&self) -> Matrix4<f32>{
        let res: Matrix4<f32> = Matrix4::from(Matrix3::from(self.rotation)); // rot_mat

        // res * scale_mat
        let scale_mat = Matrix4::from(self.scale_factor);

        // res.m[0][0] *= self.scale_factor.0;
        // res.m[1][0] *= self.scale_factor.0;
        // res.m[2][0] *= self.scale_factor.0;

        // res.m[0][1] *= self.scale_factor.1;
        // res.m[1][1] *= self.scale_factor.1;
        // res.m[2][1] *= self.scale_factor.1;
        
        // res.m[0][2] *= self.scale_factor.2;
        // res.m[1][2] *= self.scale_factor.2;
        // res.m[2][2] *= self.scale_factor.2;
        
        let translation_mat = Matrix4::new(
            [1.0, 0.0, 0.0, self.center.0],
            [0.0, 1.0, 0.0, self.center.1],
            [0.0, 0.0, 1.0, self.center.2],
            [0.0, 0.0, 0.0, 1.0],
        );
        // // translation_mat * res
        // res.m[0][3] = self.center.0;
        // res.m[1][3] = self.center.1;
        // res.m[2][3] = self.center.2;
        
        translation_mat * res * scale_mat
    }

    /// translate by a Vector3
    pub fn translate(&mut self, vector: Vector3<f32>){
        self.center = self.center + vector;
    }

    /// rotate by a quaternion
    pub fn rotate(&mut self, quaternion: Quaternion<f32>){ // where the vector3 is normalized is quaternion normalized
        self.rotation = quaternion.unit_quaternion() * self.rotation.unit_quaternion();
    }

    /// offset = camera offset from center of space
    pub fn look_at(&mut self, look_at: Vector3<f32>){
        self.rotation = Quaternion::camera_look_at_xy(self.center, look_at);
    }
}

impl Add for Space {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.join(rhs)
    }
}

impl Sub for Space {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.join_reverse(rhs)
    }
}