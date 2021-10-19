use core::fmt;
use std::ops::{Add, Div, Mul, Sub};

use crate::Zero;

use super::with_var::WithVar;
use super::one_unknown::Equa1;


#[derive(Clone, Copy)]
pub(crate) struct Equa2<T>(pub T, pub WithVar<T>, pub WithVar<T>);

impl<T> fmt::Debug for Equa2<T>
where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} + {:?} + {:?}", self.0, self.1, self.2)
    }
}

impl<T> Equa2<T>{
    pub fn new(non_var: T, a: WithVar<T>, b: WithVar<T>) -> Self {
        Equa2(non_var, a, b)
    }

    pub fn set_eq(&self, other: Self) -> ((char, T), (char, T))
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + fmt::Debug + Zero + Copy { // vars need to line up not great
        println!("{:?} = 0.0\n{:?} = 0.0", self, other);
        let var_1: (char, T) = (self.1.var, self.repl_var( Equa1::new(other.0 - self.0, other.1 - self.1) / (self.2 - other.2).val, self.2.var).solve(T::ZERO));
        let var_2: (char, T) = (self.2.var, self.plug_in(var_1.0, var_1.1).solve(T::ZERO));
        println!("{:?} = {:?}\n{:?} = {:?}", var_1.0, var_1.1, var_2.0, var_2.1);
        (var_1, var_2)
    }

    pub fn repl_var(&self, equa_1: Equa1<T>, var: char) -> Equa1<T>
    where T: Mul<T, Output = T> + Add<T, Output = T> + Copy {
        let (val_mul_of_var, other_mul_1) = self.get_var_mul_and_other(var); // 0 has no var not included
        let equa_1 = equa_1 * val_mul_of_var.val;
        Equa1::new(self.0 + equa_1.0, other_mul_1 + equa_1.1)
    }

    pub fn get_var_mul_and_other(&self, var: char) -> (WithVar<T>, WithVar<T>)
    where T: Copy{
        if self.1.var == var {
            (self.1, self.2)
        }else {
            (self.2, self.1)
        }
    }

    pub fn plug_in(&self, var: char, val: T) -> Equa1<T>
    where T: Mul<T, Output = T> + Add<T, Output = T> + Copy {
        if self.1.var == var {
            Equa1::new(self.0 + self.1.val * val, self.2)
        }else {
            Equa1::new(self.0 + self.2.val * val, self.1)
        }
    }
}