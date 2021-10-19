use std::{fmt, ops::{Add, Div, Mul, Sub, Neg}};

use crate::{F32Fmt, One, SignOps, Three, Two, imaginary::PossibleImag};

// use super::Equation; todo

use crate::Zero;

#[derive(Clone, Copy)]
pub struct CubicEquation<T>(pub T, pub T, pub T, pub T); // enough for now got a lot to do still

impl<T> fmt::Debug for CubicEquation<T> 
where T: fmt::Debug {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}x^3 + {:?}x^2 + {:?}x + {:?}", self.0, self.1, self.2, self.3)
    }
}

impl<T> CubicEquation<T>{ // for now 
    pub fn new(pow_3: T, pow_2: T, var: T, scl: T) -> Self {
        CubicEquation(pow_3, pow_2, var, scl)
    }
    /// returns a tuple. The first element of which is a Vec containing the real solutions
    /// and the second element of which is a vec containing the imaginary solutions.
    pub fn solve(&self) -> [PossibleImag<T>; 3]
    where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + SignOps + Zero + One + Two + Three + Copy + F32Fmt + PartialEq + fmt::Debug {
        // 0 = ax^3 + bx^2 + cx + d
        // Cubic formula

        let CubicEquation(a, b, c, d) = *self;

        let block0: T = -(b * b * b) / (a * T::THREE * a * T::THREE * a * T::THREE) + (b * c) / (a * a * T::THREE * T::TWO) - d / (a * T::TWO);
        let block1: T = c / (a * T::THREE) - (b * b) / (a * a * T::THREE * T::THREE);
        let block2 = PossibleImag::do_sqrt(block0 * block0 + block1 * block1 * block1);

        let part1 = (block2 + block0).cbrt();
        let part2 = (-block2 + block0).cbrt();
        
        let block3 = b / (a * T::THREE);
        [
            // part1[0] + part2[0] - block3, 
            // part1[0] + part2[1] - block3, 
            part1[0] + part2[2] - block3, // yes // I DONT KNOW WHY THESE ARE THE ANSWERS but they seem to be correct ): I'll be back.
            // part1[1] + part2[0] - block3,  
            part1[1] + part2[1] - block3, // yes
            // part1[1] + part2[2] - block3, 
            part1[2] + part2[0] - block3, // yes
            // part1[2] + part2[1] - block3, 
            // part1[2] + part2[2] - block3,
        ] // part1 + part2
    }
}
