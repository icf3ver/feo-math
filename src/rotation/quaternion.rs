//! Construct that represents a quaternion.
//! 
//! TODO
//! 
use std::{fmt, ops::{Add, Div, Mul, Neg, Rem, Sub}};

use crate::{Construct, F32Fmt, One, /* Product ,*/ SignOps, Two, Zero, axes::{self, Axes}, imaginary::ImaginaryConstruct, linear_algebra::{
        vector4::Vector4,
        vector3::Vector3,
        matrix3::Matrix3
    }};

use super::{Rotation, euler::Euler, rotor::Rotor};

// Rules
// 1   i   j   k (second)
// i   -1  k   -j
// j   -k  -1  i  
// k    j  -i  -1
// (first) ex: ij = k

/// A Quaternion with the format q = (q.0)i + (q.1)j + (q.2)k + (q.3)
/// Quaternions are usually stored 
/// `q = r + ai + bj + ck`
/// in this case however they are stored in the format 
/// `q = ai + bj + ck + r`
/// to emphasize the fact that you can essentially use the r value 
/// as a scaler.
#[derive(Clone, Copy)]
pub struct Quaternion<T>(pub T, pub T, pub T, pub T);

impl<T> PartialEq for Quaternion<T> 
where T: Neg<Output = T> + PartialEq + Copy {
    fn eq(&self, other: &Self) -> bool {
        let neg = -*other;
        (self.0 == other.0 && self.1 == other.1 && self.2 == other.2 && self.3 == other.3) ||
        (self.0 == neg.0 && self.1 == neg.1 && self.2 == neg.2 && self.3 == neg.3) 
    }
}

impl<T> Quaternion<T> {
    /// Creates a new Quaternion with defined a, b, c, and r values.
    /// # Attributes
    /// * `a` - value multiplied by the imaginary number `i` in `q = ai + bj + ck + r`
    /// * `b` - value multiplied by the imaginary number `j` in `q = ai + bj + ck + r`
    /// * `c` - value multiplied by the imaginary number `k` in `q = ai + bj + ck + r`
    /// * `r` - The real number denoted `r` in `q = ai + bj + ck + r`
    pub fn new(a: T, b: T, c: T, r: T) -> Self {
        Quaternion(a, b, c, r)
    }

    /// makes a quaternion from a Vector3 and a given real `r` value
    /// # Arguments
    /// * `vector` - A `Vector3` which maps out to (a, b, c)
    /// * `real` - The r value of the Quaternion
    pub fn new_vector_real(vector: Vector3<T>, real: T) -> Self {
        Quaternion(vector.0, vector.1, vector.2, real)
    }

    /// creates a new quaternion that corresponds to a given rotation around a given axis
    /// # Arguments
    /// * `axis` - The `Vector3` rotation axis
    /// * `angle` - the angle to rotate around that axis
    pub fn new_axis_angle(axis: Vector3<T>, angle: T) -> Self
    where T: Mul<T, Output = T> + Div<T, Output = T> + Add<T, Output = T> + F32Fmt + fmt::Debug + One + Two + Copy { // angle is in radians and axis is normalized
        let normalized_axis = axis.normalize(None);
        //println!("{:?} {:?}", normalized_axis, angle);
        let half_angle = angle / T::TWO;
        Quaternion(
            T::sin_mul(half_angle, normalized_axis.0),
            T::sin_mul(half_angle, normalized_axis.1),
            T::sin_mul(half_angle, normalized_axis.2),
            T::cos_mul(half_angle, T::ONE)
        )
    }

    /// Returns the identity quaternion that is used for multiplying 
    /// # Arguments
    /// * `zero` - the T equivalent of zero
    /// * `one` - the T equivalent of one
    pub fn identity() -> Self
    where T: Zero + One {
        Quaternion::ONE // multiplication
    }

    /// Returns the identity quaternion that is used for adding
    /// # Arguments
    /// * `zero` - the T equivalent of zero
    pub fn identity_add() -> Self
    where T: Zero {
        Quaternion::ZERO // addition
    }

    /// Returns the conjugate of the quaternion denoted `q*`
    /// # Arguments
    /// # Examples
    /// ```ignore
    /// todo!();
    /// ```
    pub fn conjugate(&self) -> Self
    where T: Neg<Output = T> + Copy {
        Quaternion(-self.0, -self.1, -self.2, self.3)
    }

    /// Returns the magnitude of the quaternion denoted `||q||`
    /// # Arguments
    /// # Examples
    /// ```ignore
    /// todo!();
    /// ```
    pub fn norm(&self) -> T
    where T: Add<T, Output = T> + Mul<T, Output = T> + F32Fmt + Copy {
        //println!("{} {}", (self.0 * self.0 + self.1 * self.1 + self.2 * self.2 + self.3 * self.3).into(), (self.0 * self.0 + self.1 * self.1 + self.2 * self.2 + self.3 * self.3).into().sqrt());
        // conjugate
        (self.0 * self.0 + self.1 * self.1 + self.2 * self.2 + self.3 * self.3).sqrt()
    }

    /// Returns the unit quaternion form of the quaternion denoted `Uq`
    /// # Arguments
    /// # Examples
    /// ```ignore
    /// todo!();
    /// ```
    pub fn unit_quaternion(&self) -> Self
    where T: Add<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + F32Fmt + Copy {
        *self / self.norm()
    }

    /// Returns the reciprocal of a quaternion denoted `q⁻¹`
    /// # Arguments
    /// # Examples
    /// ```ignore
    /// todo!();
    /// ```
    pub fn reciprocal(&self) -> Self
    where T: Add<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Copy + One {
        self.conjugate() / (self.0 * self.0 + self.1 * self.1 + self.2 * self.2 + self.3 * self.3)
    }
    
    pub fn slerp(self, to: Self, time_stamp: T::F32Fmt) -> Self 
    where T: F32Fmt + One + Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Neg<Output = T> + Copy + Div<T, Output = T>{
        let self_conj = self.conjugate();
        let cos_half_angle = (self_conj.3 * to.3 - self_conj.0 * to.0 - self_conj.1 * to.1 - self_conj.2 * to.2).intoF32Fmt();
        if cos_half_angle.abs() == T::F32Fmt::ONE {
            return self;
        }
        let half_angle = cos_half_angle.acos_mul(T::F32Fmt::ONE);

        Self::fromF32Fmt(
            (
                self.intoF32Fmt() * ((T::F32Fmt::ONE - time_stamp) * half_angle).sin_mul(T::F32Fmt::ONE)
                + to.intoF32Fmt() * (time_stamp * half_angle)
            ) / (T::F32Fmt::ONE - cos_half_angle).sqrt()
        )
    }

    pub fn sin(self) -> Self
    where T: Copy + F32Fmt + One + Sub<T, Output = T> + Add<T, Output = T> + Mul<T, Output = T> + Neg<Output = T>{
        // sin(a + bi + cj + dk)
        // sin(a + bi)cos(cj + dk) + cos(a + bi)sin(cj + dk)
        // (sin(a)cos(bi) + cos(a)sin(bi))(cos(cj)cos(dk) - sin(cj)sin(dk)) + (cos(a)cos(bi) - sin(a)sin(bi))(sin(cj)cos(dk) + cos(cj)sin(dk))
        // (sin(a)cosh(b) + icos(a)sinh(b))(cosh(c)cosh(d) - jksinh(c)sinh(d)) + (cos(a)cosh(b) - isin(a)sinh(b))(jsinh(c)cosh(d) + kcosh(c)sinh(d))
        // (sin(a)cosh(b) + icos(a)sinh(b))(cosh(c)cosh(d) - isinh(c)sinh(d)) + (cos(a)cosh(b) - isin(a)sinh(b))(jsinh(c)cosh(d) + kcosh(c)sinh(d))
        //
        // (sin(a)cosh(b)cosh(c)cosh(d) - sin(a)cosh(b)isinh(c)sinh(d) + icos(a)sinh(b)cosh(c)cosh(d) - icos(a)sinh(b)isinh(c)sinh(d))
        // + (cos(a)cosh(b)jsinh(c)cosh(d) + cos(a)cosh(b)kcosh(c)sinh(d) - isin(a)sinh(b)jsinh(c)cosh(d) - isin(a)sinh(b)kcosh(c)sinh(d))
        //
        // (sin(a)cosh(b)cosh(c)cosh(d) - isin(a)cosh(b)sinh(c)sinh(d) + icos(a)sinh(b)cosh(c)cosh(d) - iicos(a)sinh(b)sinh(c)sinh(d))
        // + (jcos(a)cosh(b)sinh(c)cosh(d) + kcos(a)cosh(b)cosh(c)sinh(d) - ijsin(a)sinh(b)sinh(c)cosh(d) - iksin(a)sinh(b)cosh(c)sinh(d))
        //
        // (sin(a)cosh(b)cosh(c)cosh(d) - isin(a)cosh(b)sinh(c)sinh(d) + icos(a)sinh(b)cosh(c)cosh(d) + cos(a)sinh(b)sinh(c)sinh(d))
        // + (jcos(a)cosh(b)sinh(c)cosh(d) + kcos(a)cosh(b)cosh(c)sinh(d) - ksin(a)sinh(b)sinh(c)cosh(d) + jsin(a)sinh(b)cosh(c)sinh(d))
        //
        // ((sin(a)cosh(b)cosh(c)cosh(d) + cos(a)sinh(b)sinh(c)sinh(d)) - i(sin(a)cosh(b)sinh(c)sinh(d) + cos(a)sinh(b)cosh(c)cosh(d)))
        // + (j(cos(a)cosh(b)sinh(c)cosh(d) + sin(a)sinh(b)cosh(c)sinh(d)) + k(cos(a)cosh(b)cosh(c)sinh(d) - sin(a)sinh(b)sinh(c)cosh(d)))
        //
        // (sin(a)cosh(b)cosh(c)cosh(d) + cos(a)sinh(b)sinh(c)sinh(d)) - i(sin(a)cosh(b)sinh(c)sinh(d) + cos(a)sinh(b)cosh(c)cosh(d))
        // + j(cos(a)cosh(b)sinh(c)cosh(d) + sin(a)sinh(b)cosh(c)sinh(d)) + k(cos(a)cosh(b)cosh(c)sinh(d) - sin(a)sinh(b)sinh(c)cosh(d))
        //
        // (sin(a)cosh(b)cosh(c)cosh(d) + cos(a)sinh(b)sinh(c)sinh(d))
        // - i(sin(a)cosh(b)sinh(c)sinh(d) + cos(a)sinh(b)cosh(c)cosh(d))
        // + j(cos(a)cosh(b)sinh(c)cosh(d) + sin(a)sinh(b)cosh(c)sinh(d))
        // + k(cos(a)cosh(b)cosh(c)sinh(d) - sin(a)sinh(b)sinh(c)cosh(d))
        let a = self.3; // r
        let b = self.0;
        let c = self.1;
        let d = self.2;
        Quaternion::new(
            a.sin_mul(b.cosh_mul(c.cosh_mul(d.cosh_mul(T::ONE)))) + a.cos_mul(b.sinh_mul(c.sinh_mul(d.sinh_mul(T::ONE)))),
            -a.sin_mul(b.cosh_mul(c.sinh_mul(d.sinh_mul(T::ONE)))) + a.cos_mul(b.sinh_mul(c.cosh_mul(d.cosh_mul(T::ONE)))),
            a.cos_mul(b.cosh_mul(c.sinh_mul(d.cosh_mul(T::ONE)))) + a.sin_mul(b.sinh_mul(c.cosh_mul(d.sinh_mul(T::ONE)))),
            a.cos_mul(b.cosh_mul(c.cosh_mul(d.sinh_mul(T::ONE)))) - a.sin_mul(b.sinh_mul(c.sinh_mul(d.cosh_mul(T::ONE)))),
        )
    }

    pub fn cos(self) -> Self{
        //
        todo!();
    }

    pub fn tan(self) -> Self{
        todo!();
    }
} 

impl Quaternion<f32> {
    /// testing
    pub fn camera_look_at_v1(pos: Vector3<f32>, look_at: Vector3<f32>) -> Self { // not done
        let frd_camrel = (look_at - pos).normalize(None); 
        let right_camrel = Vector3::cross_product(axes::NORMAL_Y_AXIS, frd_camrel).normalize(None); // axis 2
        
        Quaternion::from_axes(right_camrel, frd_camrel)
    }

    /// testing
    pub fn from_axes(left: Vector3<f32>, frd: Vector3<f32>) -> Self { // merge 
        let up = Vector3::cross_product( frd, left); //.normalize(None); // axis 2

        // so the rotation matrix is:
        // [  left ] this is row representation!
        // [  up   ] each axis is a row in the 
        // [  frd  ] corresponding rotation matrix

        let derived_rot_mat = Matrix3::new( // rmb in col format
            left.into(),
            up.into(),
            frd.into(),
        ).transpose();
        
        derived_rot_mat.into()
    }
    
    /// Camera Objects look along the -Z axis so the look_at function for a camera object
    /// needs to be modified a little. flipping the pos and the look_at should do the trick.
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
    /// # Examples
    /// ```ignore //rust
    /// use crate::linear_algebra::Vector3;
    /// let pos = Vector3::new(0.0, 0.0, 0.0);
    /// let look_at = Vector3::new(1.0, 1.0, 1.0);
    /// assert_eq!(camera_look_at(pos, look_at), todo!());
    /// ```
    pub fn camera_look_at(pos: Vector3<f32>, look_at: Vector3<f32>) -> Quaternion<f32> {
        Quaternion::look_at_xy(look_at, pos)
    }

    pub fn f32_identity() -> Self{
        Quaternion::new( 0.0_f32, 0.0_f32, 0.0_f32, 1.0_f32)
    }
}

impl<T> Construct<T> for Quaternion<T> where T: Construct<T> {}
impl<T> ImaginaryConstruct<T> for Quaternion<T> where T: Construct<T> {}
impl<T> Rotation<T> for Quaternion<T> where T: Construct<T> {
    /// # Examples
    /// ```ignore //rust
    /// use crate::linear_algebra::Vector3;
    /// let pos = Vector3::new(0.0, 0.0, 0.0);
    /// let look_at = Vector3::new(1.0, 1.0, 1.0);
    /// assert_eq!(camera_look_at(pos, look_at), Quaternion::new());
    /// ```
    fn look_at_xy(pos: Vector3<T>, look_at: Vector3<T>) -> Self {
        let frd_camrel = (look_at - pos).normalize(None);
        let right_camrel = Vector3::cross_product(Self::Y_AXIS, frd_camrel).normalize(None);

        let frd_no_y = Vector3::cross_product(right_camrel, Self::Y_AXIS).normalize(None);//. z;

        //             0 when in line with each other
        //            -y axis 
        // -sqrt(2)/2 _|_ sqrt(2)/2
        //          ___|___ frd 1
        // -sqrt(2)/2 _|_ sqrt(2)/2
        //             |
        //             0 when exactly opposite
        let self_x_axis_turn = T::asin_mul(
            Vector3::dot_product(-Self::Y_AXIS, frd_camrel),
            T::ONE
        );

        let self_y_axis_turn = T::atan2_mul(
            Vector3::dot_product(Self::X_AXIS, frd_no_y),
            Vector3::dot_product(Self::X_AXIS, right_camrel), 
            T::ONE
        );

        Quaternion::new_axis_angle(right_camrel, self_x_axis_turn) *
        Quaternion::new_axis_angle(Self::Y_AXIS, self_y_axis_turn)
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

    /// # Examples
    /// ```ignore //rust
    /// use crate::linear_algebra::Vector3;
    /// let pos = Vector3::new(0.0, 0.0, 0.0);
    /// let look_at = Vector3::new(1.0, 1.0, 1.0);
    /// assert_eq!(camera_look_at(pos, look_at), todo!());
    /// ```
    fn camera_look_at_xy(pos: Vector3<T>, look_at: Vector3<T>) -> Self { Quaternion::look_at_xy(look_at, pos) }

    fn camera_look_at_xz(pos: Vector3<T>, look_at: Vector3<T>) -> Self { Quaternion::look_at_xz(look_at, pos) }

    fn camera_look_at_yz(pos: Vector3<T>, look_at: Vector3<T>) -> Self { Quaternion::look_at_yz(look_at, pos) }

    fn camera_look_at_lock(pos: Vector3<T>, look_at: Vector3<T>, locked_axis: Vector3<T>) -> Self { 
        Quaternion::look_at_lock(look_at, pos, locked_axis)
    }
}
impl<T> F32Fmt for Quaternion<T> 
where T: F32Fmt + Copy + Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + One { 
    type F32Fmt = Quaternion<T::F32Fmt>;
    #[inline]
    fn intoF32Fmt(self) -> Self::F32Fmt {
        Quaternion(self.0.intoF32Fmt(), self.1.intoF32Fmt(), self.2.intoF32Fmt(), self.3.intoF32Fmt())
    } 
    #[inline]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self {
        let vec = &f32_fmt;
        Quaternion(T::fromF32Fmt(vec.0), T::fromF32Fmt(vec.1), T::fromF32Fmt(f32_fmt.2), T::fromF32Fmt(f32_fmt.3))
    }
 
    fn sqrt(self) -> Self {
        let v_part = Vector3(self.0, self.1, self.2);
        let r_part = self.3;

        let block0 = (self.norm() + r_part).f32_const_mul(0.5).sqrt();

        Quaternion::new_vector_real(v_part.unit_vector() * block0, block0)
    }

    fn cbrt(self) -> Self {
        todo!();
    }

    fn f32_const_mul(self, constant: f32) -> Self {
        Quaternion(
            self.0.f32_const_mul(constant), 
            self.1.f32_const_mul(constant), 
            self.2.f32_const_mul(constant), 
            self.3.f32_const_mul(constant)
        )
    }

    fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        self.sin() * mul_by
    }

    fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        self.cos() * mul_by
    }

    fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized {
        self.tan() * mul_by
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

impl<T> From<Matrix3<T>> for Quaternion<T> 
where T: F32Fmt + Add<T, Output = T> + Sub<T, Output = T> + Div<T, Output = T> + Mul<T, Output = T> + fmt::Debug + Two + One + Zero + Copy {
    fn from(mat3: Matrix3<T>) -> Self {
        // get axis
        let axis = mat3.eigen_vector(T::ONE); // nd

        // get angle
        // |\ 2
        // |_\
        //  1
        
        let trace = mat3.trace();
        let angle = T::acos_mul((trace - T::ONE) / T::TWO, T::ONE);

        Quaternion::new_axis_angle(axis, angle)
    }
}
impl<T> From<Axes<T>> for Quaternion<T> 
where T: Construct<T> + Copy {
    fn from(a: Axes<T>) -> Self {
        Into::<Matrix3<T>>::into(a).into()
    }
}
impl<T> From<Euler<T>> for Quaternion<T> 
where T: Construct<T> + Copy {
    fn from(e: Euler<T>) -> Self {
        Into::<Matrix3<T>>::into(e).into()
    }
}
impl<T> From<Rotor<T>> for Quaternion<T> 
where T: Construct<T> + Copy {
    fn from(r: Rotor<T>) -> Self {
        Into::<Matrix3<T>>::into(r).into()
    }
}

impl<T> From<Quaternion<T>> for Matrix3<T> 
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Copy + One + Two {
    fn from(other: Quaternion<T>) -> Matrix3<T> {
        let i = other.0; let i_squared = i * i;
        let j = other.1; let j_squared = j * j;
        let k = other.2; let k_squared = k * k;
        let r = other.3; // distance from center 

        // p'
        
        let two: T = T::TWO;
        let one: T = T::ONE;
        // let s: T = one / (i_squared + j_squared + k_squared + r * r);
        // Matrix3::new(
        //     [one - two * s * (j_squared + k_squared),               two * s * (i * j - k * r),               two * s * (i * k + j * r)],
        //     [              two * s * (i * j + k * r), one - two * s * (i_squared + k_squared),               two * s * (j * k - i * r)],
        //     [              two * s * (i * k - j * r),               two * s * (j * k + i * r), one - two * s * (i_squared + j_squared)],
        // )

        // let s: T = two / (i_squared + j_squared + k_squared + r * r);
        // Matrix3::new(
        //     [one - s * (j_squared + k_squared),               s * (i * j - k * r),               s * (i * k + j * r)],
        //     [              s * (i * j + k * r), one - s * (i_squared + k_squared),               s * (j * k - i * r)],
        //     [              s * (i * k - j * r),               s * (j * k + i * r), one - s * (i_squared + j_squared)],
        // )


        let two_s: T = two / (i_squared + j_squared + k_squared + r * r);

        Matrix3::new(
            [one - two_s * (j_squared + k_squared),               two_s * (i * j - k * r),               two_s * (i * k + j * r)],
            [              two_s * (i * j + k * r), one - two_s * (i_squared + k_squared),               two_s * (j * k - i * r)],
            [              two_s * (i * k - j * r),               two_s * (j * k + i * r), one - two_s * (i_squared + j_squared)],
        )
    }
}
impl<T> From<Quaternion<T>> for Vector4<T> where T: Copy { fn from(other: Quaternion<T>) -> Vector4<T> { Vector4(other.0, other.1, other.2, other.3) } } // rmv pointless
impl<T> From<Quaternion<T>> for Vector3<T> where T: Copy { fn from(other: Quaternion<T>) -> Vector3<T> { Vector3(other.0, other.1, other.2) } } // just the vector/imaginary part
impl<T> From<Quaternion<T>> for [T; 4] where T: Copy { fn from(other: Quaternion<T>) -> [T; 4] { [other.0, other.1, other.2, other.3] } }
impl<T> fmt::Debug for Quaternion<T> where T: Copy + fmt::Debug {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}i + {:?}j + {:?}k + {:?}", self.0, self.1, self.2, self.3)
    }
}

impl<T> Zero for Quaternion<T> where T: Zero { const ZERO: Self = Quaternion(T::ZERO, T::ZERO, T::ZERO, T::ZERO); }
impl<T> One for Quaternion<T> where T: One + Zero { const ONE: Self = Quaternion(T::ZERO, T::ZERO, T::ZERO, T::ONE); }
impl<T> Two for Quaternion<T> where T: Two + Zero { const TWO: Self = Quaternion(T::ZERO, T::ZERO, T::ZERO, T::TWO); }


impl<T> Add for Quaternion<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Quaternion(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2, self.3 + rhs.3)
    }
}
impl<T> Add<T> for Quaternion<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        Quaternion(self.0 + rhs, self.1 + rhs, self.2 + rhs, self.3 + rhs)
    }
}

impl<T> Sub for Quaternion<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Quaternion(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2, self.3 - rhs.3)
    }
}
impl<T> Sub<T> for Quaternion<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        Quaternion(self.0 - rhs, self.1 - rhs, self.2 - rhs, self.3 - rhs)
    }
}

impl<T> Mul for Quaternion<T> 
where T: Sub<T, Output = T> + Mul<T, Output = T> + Add<T, Output = T> + Copy {
    type Output = Self;

    // hamilton product
    fn mul(self, rhs: Self) -> Self::Output {
        // rhs = rhs.3 + rhs.0i+ rhs.1j + rhs.2k
        // self = self.3 + self.0i + self.1j + self.2k
        //
        // Rules
        // 1   i   j   k (second)
        // i   -1  k   -j
        // j   -k  -1  i  
        // k    j  -i  -1
        // (first) ex: ij = k
        //
        
        let r = self.3 * rhs.3 - self.0 * rhs.0 - self.1 * rhs.1 - self.2 * rhs.2;
        let i = self.3 * rhs.0 + self.0 * rhs.3 + self.1 * rhs.2 - self.2 * rhs.1;
        let j = self.3 * rhs.1 - self.0 * rhs.2 + self.1 * rhs.3 + self.2 * rhs.0;
        let k = self.3 * rhs.2 + self.0 * rhs.1 - self.1 * rhs.0 + self.2 * rhs.3;

        Quaternion(i, j, k, r)
        
        // let v1: Vector3<T> = self.into();
        // let v2: Vector3<T> = rhs.into();

        // let r1 = self.3;
        // let r2 = rhs.3;

        // let vector_part: Vector3<T> = v2 * r1 + v1 * r2 + Vector3::cross_product(v1, v2);
        // Quaternion(vector_part.0, vector_part.1, vector_part.2, r1 * r2 - Vector3::dot_product(v1, v2))
    }
}
impl<T> Mul<T> for Quaternion<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Quaternion(self.0 * rhs, self.1 * rhs, self.2 * rhs, self.3 * rhs)
    }
}

impl<T> Div for Quaternion<T>
where T: Copy + Neg<Output = T> + Sub<T, Output = T> + Add<T, Output = T> + Div<T, Output = T> {
    type Output = Self;

    #[allow(clippy::many_single_char_names)]
    fn div(self, rhs: Self) -> Self::Output {
        // q1 = a + bi + cj + dk
        // q2 = e + fi + gj + hk
        
        // q1/q2 = a/e + a/fi + a/gj + a/hk
        //       + bi/e + bi/fi + bi/gj + bi/hk
        //       + cj/e + cj/fi + cj/gj + cj/hk
        //       + dk/e + dk/fi + dk/gj + dk/hk
        
        // q1/q2 = a/e + a/fi + a/gj + a/hk
        //       + bi/e + b/f + bi/gj + bi/hk
        //       + cj/e + cj/fi + c/g + cj/hk
        //       + dk/e + dk/fi + dk/gj + d/h
        
        // q1/q2 = a/e + ai/fii + aj/gjj + ak/hkk
        //       + bi/e + b/f + bij/gjj + bik/hkk
        //       + cj/e + cji/fii + c/g + cjk/hkk
        //       + dk/e + dki/fii + dkj/gjj + d/h
        
        // q1/q2 = a/e - ai/f - aj/g - ak/h
        //       + bi/e + b/f - bij/g - bik/h
        //       + cj/e - cji/f + c/g - cjk/h
        //       + dk/e - dki/f - dkj/g + d/h
        
        // q1/q2 = a/e - ai/f - aj/g - ak/h
        //       + bi/e + b/f - bk/g + bj/h
        //       + cj/e + ck/f + c/g - ci/h
        //       + dk/e - dj/f + di/g + d/h
        
        // q1/q2 = a/e - (a/f)i - (a/g)j - (a/h)k
        //       + (b/e)i + b/f - (b/g)k + (b/h)j
        //       + (c/e)j + (c/f)k + c/g - (c/h)i
        //       + (d/e)k - (d/f)j + (d/g)i + d/h
        
        // q1/q2 = (a/e + b/f + c/g + d/h) 
        //       + (-a/f + b/e - c/h + d/g)i
        //       + (-a/g + b/h + c/e - d/f)j 
        //       + (-a/h - b/g + c/f + d/e)k
        
        let Quaternion(a, b, c, d) = self;
        let Quaternion(e, f, g, h) = rhs;

        Quaternion(
            -a/f + b/e - c/h + d/g,
            -a/g + b/h + c/e - d/f,
            -a/h - b/g + c/f + d/e,
            a/e + b/f + c/g + d/h
        )
    }
}
impl<T> Div<T> for Quaternion<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Quaternion(self.0 / rhs, self.1 / rhs, self.2 / rhs, self.3 / rhs)
    }
}

impl<T> Rem for Quaternion<T> 
where T: Copy + Neg<Output = T> + Sub<T, Output = T> + Add<T, Output = T> + Rem<T, Output = T> {
    type Output = Self;

    /// TODO Test
    #[allow(clippy::many_single_char_names)]
    fn rem(self, rhs: Self) -> Self::Output {
        // same as div

        let Quaternion(a, b, c, d) = self;
        let Quaternion(e, f, g, h) = rhs;

        Quaternion(
            -a%f + b%e - c%h + d%g,
            -a%g + b%h + c%e - d%f,
            -a%h - b%g + c%f + d%e,
            a%e + b%f + c%g + d%h
        )
    }
}

impl<T> Rem<T> for Quaternion<T> 
where T: Rem<T, Output = T> + Copy {
    type Output = Self;

    fn rem(self, rhs: T) -> Self::Output {
        Quaternion(self.0 % rhs, self.1 % rhs, self.2 % rhs, self.3 % rhs)
    }
}

impl<T> Neg for Quaternion<T> 
where T: Neg<Output = T> + Copy {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Quaternion(-self.0, -self.1, -self.2, -self.3)
    }
}

impl<T> SignOps for Quaternion<T> 
where T: SignOps + Add<T, Output = T> {
    fn ptcopysign(self, sign: Self) -> Self {
        Quaternion(
            self.0.ptcopysign(sign.0), 
            self.1.ptcopysign(sign.1), 
            self.2.ptcopysign(sign.2), 
            self.3.ptcopysign(sign.3)
        )
    }
    
    // Does not really make sense
    fn ptsignum(self) -> i8 {
        (self.0 + self.1 + self.2 + self.3).ptsignum()
    }

    fn abs(self) -> Self {
        Quaternion(
            self.0.abs(), 
            self.1.abs(), 
            self.2.abs(), 
            self.3.abs()
        )
    }
}