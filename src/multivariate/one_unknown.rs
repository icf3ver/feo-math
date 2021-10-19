use std::{ ops::{Div, Mul, Sub}};

use super::with_var::WithVar;

#[derive(Clone, Copy)]
pub(crate) struct Equa1<T>(pub T, pub WithVar<T>);

impl<T> Equa1<T>{
    pub fn new(non_var: T, a: WithVar<T>) -> Self {
        Equa1(non_var, a)
    }
    #[allow(dead_code)]
    pub fn set_eq(&self, other: Self) -> T 
    where T: Div<T, Output = T> + Sub<T, Output = T> + Copy {
        (self.0 - other.0) / (other.1 - self.1).val
    }
    pub fn solve(&self, other: T) -> T
    where T: Sub<T, Output = T> + Div<T, Output = T> + Copy {
        (other - self.0) / (self.1).val
    }
    pub fn get_var(&self, var: char) -> Option<WithVar<T>>
    where T: Copy {
        if self.1.var == var {
            return Some(self.1);
        }
        println!("{} {}", self.1.var, var);
        None
    }
}

impl<T> Mul<T> for Equa1<T>
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Equa1(self.0 * rhs, self.1 * rhs)
    }
}

impl<T> Mul<WithVar<T>> for Equa1<T>
where T: Mul<T, Output = T> + Copy {
    type Output = Result<Self, ()>;

    fn mul(self, rhs: WithVar<T>) -> Self::Output {
        match self.get_var(rhs.var){ // not good 
            Some(_var) => Ok(Equa1(self.0, self.1 * rhs)),
            None => Err(())
        }
    }
}

impl<T> Div<WithVar<T>> for Equa1<T>
where T: Div<T, Output = T> + Copy {
    type Output = Result<Self, ()>;

    fn div(self, rhs: WithVar<T>) -> Self::Output {
        match self.get_var(rhs.var){ // not good 
            Some(_var) => Ok(Equa1(self.0, self.1 / rhs)),
            None => Err(())
        }
    }
}

impl<T> Div<T> for Equa1<T>
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Equa1(self.0/rhs, self.1 / rhs)
    }
}