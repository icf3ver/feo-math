use {
    crate::Construct,
    std::ops::Mul
};

pub mod matrix3;
pub mod matrix4;
pub mod vector3;
pub mod vector4;

pub trait Vector<T>: Construct<T> 
{ /* none */ }

pub trait SqMatrix<T, U>: Construct<T> + Mul<U, Output = U>
where T: Construct<T>, 
      U: Vector<T>,
{ /* none */ }