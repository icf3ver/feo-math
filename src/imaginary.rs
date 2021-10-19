use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::{F32Fmt, One, SignOps, Three, Two, Zero};

use super::Construct;

pub trait ImaginaryConstruct<T>: Construct<T> { /* todo */ }

// /// Imaginary construct type .0 n .1
// pub struct Imaginary<T>(char, T);

// pub struct Complex<T>(T, Vec<Imaginary<T>>);

/// For now
#[derive(Clone, Copy)]
pub struct Simple2DComplexNumber<T>(pub T, pub T);

impl<T> Simple2DComplexNumber<T> {
    pub fn is_real (self) -> bool 
    where T: Zero + PartialEq {
        self.1 == T::ZERO
    }

    pub fn to_real(self) -> Result<T, ()> 
    where T: Zero + PartialEq {
        if self.1 == T::ZERO {
            Ok(self.0)            
        } else {
            Err(())
        }
    }

    pub fn cbrt (self) -> [Self; 3]
    where T: F32Fmt + PartialEq + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + One + Two + Three + SignOps + Div<T, Output = T> + Copy {
        let Simple2DComplexNumber(real, imag) = self;
        let r = (real * real + imag * imag).sqrt();
        let angle_pfi = T::atan2_mul(imag, real, T::ONE);
        
        let cbrt_r = r.cbrt();
        // Answers: 
        // 1)
        // r.cbrt() * (e^(iTheta/3)) = r.cbrt() * (e^(i(Theta/3))
        // Eulers formula: (e^(i(Theta/3)) = cos(Theta/3) + sin(Theta/3)i
        // r.cbrt() * cos(Theta/3) + sin(Theta/3)i
        // r.cbrt() * cos(Theta/3) + (r.cbrt() * sin(Theta/3))i
        // 2)
        // r.cbrt() * (e^(iTheta/3 + 2iPI/3)) = r.cbrt() * (e^(i((Theta+2PI)/3))
        // ... Eulers formula
        // r.cbrt() * cos((Theta+2PI)/3) + sin((Theta+2PI)/3)i
        // (r.cbrt() * cos((Theta+2PI)/3)) + (r.cbrt() * sin((Theta+2PI)/3))i
        // 3)
        // r.cbrt() * (e^(iTheta/3 - 2iPI/3)) = r.cbrt() * (e^(i((Theta-2PI)/3))
        // ... Eulers formula
        // r.cbrt() * cos((Theta-2PI)/3) + sin((Theta-2PI)/3)i
        // (r.cbrt() * cos((Theta-2PI)/3)) + (r.cbrt() * sin((Theta-2PI)/3))i

        let block0 = angle_pfi / T::THREE;
        let answer0 = Simple2DComplexNumber(T::cos_mul(block0, cbrt_r), T::sin_mul(block0, cbrt_r));
        let block1 = (angle_pfi + T::f32_const_mul(T::TWO, std::f32::consts::PI)) / T::THREE;
        let answer1 = Simple2DComplexNumber(T::cos_mul(block1, cbrt_r), T::sin_mul(block1, cbrt_r));
        let block2 = (angle_pfi - T::f32_const_mul(T::TWO, std::f32::consts::PI)) / T::THREE;
        let answer2 = Simple2DComplexNumber(T::cos_mul(block2, cbrt_r), T::sin_mul(block2, cbrt_r));

        [answer0, answer1, answer2]
    }
}

impl<T> Add<Simple2DComplexNumber<T>> for Simple2DComplexNumber<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: Simple2DComplexNumber<T>) -> Self::Output {
        Simple2DComplexNumber(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl<T> Add<T> for Simple2DComplexNumber<T> 
where T: Add<T, Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        Simple2DComplexNumber(self.0 + rhs, self.1)
    }
}

impl<T> Sub<Simple2DComplexNumber<T>> for Simple2DComplexNumber<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: Simple2DComplexNumber<T>) -> Self::Output {
        Simple2DComplexNumber(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl<T> Sub<T> for Simple2DComplexNumber<T> 
where T: Sub<T, Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        Simple2DComplexNumber(self.0 - rhs, self.1)
    }
}

impl<T> Mul<Simple2DComplexNumber<T>> for Simple2DComplexNumber<T> 
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Sub<T, Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: Simple2DComplexNumber<T>) -> Self::Output {
        Simple2DComplexNumber(self.0 * rhs.0 - self.1 * rhs.1, self.0 * rhs.1 + self.1 * rhs.0)
    }
}

impl<T> Mul<T> for Simple2DComplexNumber<T> 
where T: Mul<T, Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Simple2DComplexNumber(self.0 * rhs, self.1 * rhs)
    }
}

impl<T> Div<Simple2DComplexNumber<T>> for Simple2DComplexNumber<T> 
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: Simple2DComplexNumber<T>) -> Self::Output {
        let numerator = (self.0 * rhs.0 +/*-*/ self.1 * /*-*/rhs.1, self.1 * rhs.0 -/*+*/ self.0 * /*-*/rhs.1);
        let real_denominator = rhs.0 * rhs.0 +/*-*/ rhs.1 * /*-*/rhs.1;
        Simple2DComplexNumber(numerator.0 / real_denominator, numerator.1 / real_denominator)
    }
}

impl<T> Div<T> for Simple2DComplexNumber<T> 
where T: Div<T, Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Simple2DComplexNumber(self.0 / rhs, self.1 / rhs)
    }
}

impl<T> Neg for Simple2DComplexNumber<T> 
where T: Neg<Output = T> + Copy {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Simple2DComplexNumber(-self.0, -self.1)
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum PossibleImag<T> {
    Complex(T, T),
    Imag(T),
    Real(T),
}
impl<T> PossibleImag<T> {
    pub fn new(real: T, imag: T) -> Self 
    where T: Zero + PartialEq {
        if imag == T::ZERO {
            PossibleImag::Real(real)
        } else if real == T::ZERO {
            PossibleImag::Imag(imag)
        } else {
            PossibleImag::Complex(real, imag)
        }
    }

    pub fn do_sqrt(other: T) -> Self 
    where T: F32Fmt + SignOps + Clone {
        let n = other.clone().abs().sqrt();
        if other.ptsignum() == -1 {
            PossibleImag::Imag(n)
        } else {
            PossibleImag::Real(n)
        }
    }

    pub fn cbrt(self) -> [Self; 3]
    where T: Zero + F32Fmt + PartialEq + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + One + Two + Three + SignOps + Div<T, Output = T> + Copy {
        let (real, imag) = match self {
            PossibleImag::Complex(r, i) => (r, i),
            PossibleImag::Imag(i) => (T::ZERO, i),
            PossibleImag::Real(r) => (r, T::ZERO)
        };
        let r = (real * real + imag * imag).sqrt();
        let angle_pfi = T::atan2_mul(imag, real, T::ONE);
        
        let cbrt_r = r.cbrt();
        // Answers: 
        // 1)
        // r.cbrt() * (e^(iTheta/3)) = r.cbrt() * (e^(i(Theta/3))
        // Eulers formula: (e^(i(Theta/3)) = cos(Theta/3) + sin(Theta/3)i
        // r.cbrt() * cos(Theta/3) + sin(Theta/3)i
        // r.cbrt() * cos(Theta/3) + (r.cbrt() * sin(Theta/3))i
        // 2)
        // r.cbrt() * (e^(iTheta/3 + 2iPI/3)) = r.cbrt() * (e^(i((Theta+2PI)/3))
        // ... Eulers formula
        // r.cbrt() * cos((Theta+2PI)/3) + sin((Theta+2PI)/3)i
        // (r.cbrt() * cos((Theta+2PI)/3)) + (r.cbrt() * sin((Theta+2PI)/3))i
        // 3)
        // r.cbrt() * (e^(iTheta/3 - 2iPI/3)) = r.cbrt() * (e^(i((Theta-2PI)/3))
        // ... Eulers formula
        // r.cbrt() * cos((Theta-2PI)/3) + sin((Theta-2PI)/3)i
        // (r.cbrt() * cos((Theta-2PI)/3)) + (r.cbrt() * sin((Theta-2PI)/3))i

        let block0 = angle_pfi/T::THREE;
        let answer0 = PossibleImag::new(T::cos_mul(block0, cbrt_r), T::sin_mul(block0, cbrt_r));
        let block1 = (angle_pfi + T::f32_const_mul(T::TWO, std::f32::consts::PI)) / T::THREE;
        let answer1 = PossibleImag::new(T::cos_mul(block1, cbrt_r), T::sin_mul(block1, cbrt_r));
        let block2 = (angle_pfi - T::f32_const_mul(T::TWO, std::f32::consts::PI)) / T::THREE;
        let answer2 = PossibleImag::new(T::cos_mul(block2, cbrt_r), T::sin_mul(block2, cbrt_r));

        [answer0, answer1, answer2]
    }
}

impl<T: Add<T, Output = T> + PartialEq + Zero> Add<PossibleImag<T>> for PossibleImag<T> {
    type Output = PossibleImag<T>;

    fn add(self, rhs: PossibleImag<T>) -> Self::Output {
        let (mut real, mut imag) = match self {
            PossibleImag::Complex(r, i) => (r, i),
            PossibleImag::Imag(i) => (T::ZERO, i),
            PossibleImag::Real(r) => (r, T::ZERO)
        };

        match rhs {
            PossibleImag::Complex(r, i) => {
                real = real + r;
                imag = imag + i;
            },
            PossibleImag::Imag(i) => imag = imag + i,
            PossibleImag::Real(r) => real = real + r
        }

        match (real, imag) {
            (real, zero) if zero == T::ZERO => PossibleImag::Real(real),
            (zero, imag) if zero == T::ZERO => PossibleImag::Imag(imag),
            (real, imag) => PossibleImag::Complex(real, imag)
        }
        
    }
}

impl<T: Add<T, Output = T>> Add<T> for PossibleImag<T> {
    type Output = PossibleImag<T>;

    fn add(self, rhs: T) -> Self::Output {
        match self {
            PossibleImag::Complex(r, i) => PossibleImag::Complex(r + rhs, i),
            PossibleImag::Imag(i) => PossibleImag::Complex(rhs, i),
            PossibleImag::Real(r) => PossibleImag::Real(r + rhs)
        }
    }
}

impl<T: Sub<T, Output = T> + PartialEq + Zero> Sub<PossibleImag<T>> for PossibleImag<T> {
    type Output = PossibleImag<T>;

    fn sub(self, rhs: PossibleImag<T>) -> Self::Output {
        let (mut real, mut imag) = match self {
            PossibleImag::Complex(r, i) => (r, i),
            PossibleImag::Imag(i) => (T::ZERO, i),
            PossibleImag::Real(r) => (r, T::ZERO)
        };

        match rhs {
            PossibleImag::Complex(r, i) => {
                real = real - r;
                imag = imag - i;
            },
            PossibleImag::Imag(i) => imag = imag - i,
            PossibleImag::Real(r) => real = real - r
        }

        match (real, imag) {
            (real, zero) if zero == T::ZERO  => PossibleImag::Real(real),
            (zero, imag) if zero == T::ZERO  => PossibleImag::Imag(imag),
            (real, imag) => PossibleImag::Complex(real, imag)
        }
    }
}

impl<T: Sub<T, Output = T> + Neg<Output = T>> Sub<T> for PossibleImag<T> {
    type Output = PossibleImag<T>;

    fn sub(self, rhs: T) -> Self::Output {
        match self {
            PossibleImag::Complex(r, i) => PossibleImag::Complex(r - rhs, i),
            PossibleImag::Imag(i) => PossibleImag::Complex(-rhs, i),
            PossibleImag::Real(r) => PossibleImag::Real(r - rhs)
        }
    }
}

impl<T> Mul<PossibleImag<T>> for PossibleImag<T> 
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Neg<Output = T> + Zero + PartialEq + Copy {
    type Output = PossibleImag<T>;

    fn mul(self, rhs: PossibleImag<T>) -> Self::Output {
        let (mut real, mut imag) = match self {
            PossibleImag::Complex(r, i) => (r, i),
            PossibleImag::Imag(i) => (T::ZERO, i),
            PossibleImag::Real(r) => (r, T::ZERO)
        };

        match rhs {
            PossibleImag::Complex(r, i) => {
                let tmp_real = real * r - imag * i;
                imag = real * i + imag * r;
                real = tmp_real;
            },
            PossibleImag::Imag(i) => {
                let tmp_real = - imag * i;
                imag = real * i;
                real = tmp_real;
            },
            PossibleImag::Real(r) => {
                real = real * r;
                imag = imag * r;
            }
        }

        match (real, imag) {
            (real, zero) if zero == T::ZERO => PossibleImag::Real(real),
            (zero, imag) if zero == T::ZERO => PossibleImag::Imag(imag),
            _ => PossibleImag::Complex(real, imag)
        }
    }
}

impl<T: Mul<T, Output = T>> Mul<T> for PossibleImag<T> 
where T: Copy {
    type Output = PossibleImag<T>;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            PossibleImag::Complex(r, i) => PossibleImag::Complex(r * rhs, i * rhs),
            PossibleImag::Imag(i) => PossibleImag::Imag(i * rhs),
            PossibleImag::Real(r) => PossibleImag::Real(r * rhs)
        }
    }
}

impl<T> Div<PossibleImag<T>> for PossibleImag<T> 
where T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T> + Neg<Output = T> + Zero + PartialEq + Copy {
    type Output = PossibleImag<T>;

    fn div(self, rhs: PossibleImag<T>) -> Self::Output {
        let (mut real, mut imag) = match self {
            PossibleImag::Complex(r, i) => (r, i),
            PossibleImag::Imag(i) => (T::ZERO, i),
            PossibleImag::Real(r) => (r, T::ZERO)
        };

        match rhs {
            PossibleImag::Complex(r, i) => {
                { // new denominator
                    let tmp_real = real * r +/*-*/ imag * /*-*/i;
                    imag = imag * r -/*+*/ real * /*-*/i;
                    real = tmp_real;
                }
                let real_denominator = r * r +/*-*/ i * /*-*/i;
                imag = imag / real_denominator;
                real = real / real_denominator;
            },
            PossibleImag::Imag(i) => {
                { // new denominator
                    let tmp_real = imag * i;
                    imag = real * i;
                    real = tmp_real;
                }
                let real_denominator = - i * i;
                imag = imag / real_denominator;
                real = real / real_denominator;
            },
            PossibleImag::Real(r) => {
                real = real / r;
                imag = imag / r;
            }
        }

        match (real, imag) {
            (real, zero) if zero == T::ZERO => PossibleImag::Real(real),
            (zero, imag) if zero == T::ZERO => PossibleImag::Imag(imag),
            _ => PossibleImag::Complex(real, imag)
        }
    }
}

impl<T: Div<T, Output = T>> Div<T> for PossibleImag<T>  
where T: Copy {
    type Output = PossibleImag<T>;

    fn div(self, rhs: T) -> Self::Output {
        match self {
            PossibleImag::Complex(r, i) => PossibleImag::Complex(r / rhs, i / rhs),
            PossibleImag::Imag(i) => PossibleImag::Imag(i / rhs),
            PossibleImag::Real(r) => PossibleImag::Real(r / rhs)
        }
    }
}

impl<T> Neg for PossibleImag<T> 
where T: Neg<Output = T> {
    type Output = PossibleImag<T>;

    fn neg(self) -> Self::Output {
        match self {
            PossibleImag::Real(r) => PossibleImag::Real(-r),
            PossibleImag::Imag(i) => PossibleImag::Imag(-i),
            PossibleImag::Complex(r, i) => PossibleImag::Complex(-r, -i),
        }
    }
}