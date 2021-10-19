use std::{fmt, ops::{Add, Div, Mul, Neg, Rem, Shl, Shr, Sub}};

pub mod multivariate;
pub mod polynomial;
pub mod imaginary;
pub mod linear_algebra;
pub mod rotation;
pub mod utils; // unfinished move forward

pub use rotation::axes;

pub trait Construct<T>: fmt::Debug + Copy + Clone + PartialEq // Is typed
    + Add<Self, Output = Self> + Add<T, Output = Self>
    + Sub<Self, Output = Self> + Sub<T, Output = Self> 
    + Mul<Self, Output = Self> + Mul<T, Output = Self> 
    + Div<Self, Output = Self> + Div<T, Output = Self>
    + Rem<Self, Output = Self> + Rem<T, Output = Self>
    + Neg<Output = Self> + SignOps
    + Zero + One + Two
    + F32Fmt {}

/// This trait allows you to set requirements for Types of typed Constructs in groups.
/// TODO explain + example + fix issue note primitive
pub trait Typed<T> {} // A concept 
impl<T: Construct<A>, A> Typed<T> for A {} // Note that when no longer typed typed is set to the last type.

// pub trait Aggregate<T>{
//     fn ag(lhs: T, rhs: T) -> T;
// }

// pub struct None(!);
// impl<T> Aggregate<T> for None {
//     fn ag(_: T, _: T) -> T {
//         panic!("Not allowed");
//     }
// }
// pub struct Product(!);
// impl<T> Aggregate<T> for Product 
// where T: Mul<T, Output=T> {
//     fn ag(lhs: T, rhs: T) -> T{
//         lhs * rhs
//     }
// }
// pub struct Sum(!);
// impl<T> Aggregate<T> for Sum 
// where T: Add<T, Output=T> {
//     fn ag(lhs: T, rhs: T) -> T {
//         lhs + rhs
//     }
// }

impl<T: fmt::Debug + Copy + Clone + PartialEq
        + Add<T, Output = T>
        + Sub<T, Output = T>
        + Mul<T, Output = T>
        + Div<T, Output = T>
        + Rem<T, Output = T>
        + Neg<Output = T> + SignOps
        + Zero + One + Two
        + F32Fmt
        + Primitive> Construct<T> for T {}

pub trait IntUtils{
    fn isqrt(self) -> Self;
    fn icbrt(self) -> Self;
    fn isin_mul(self, mul_by: Self) -> Self;
    fn icos_mul(self, mul_by: Self) -> Self;
    fn itan_mul(self, mul_by: Self) -> Self;
    fn iasin_mul(self, mul_by: Self) -> Self;
    fn iacos_mul(self, mul_by: Self) -> Self;
    fn iatan_mul(self, mul_by: Self) -> Self;
    fn iatan2_mul(self, other: Self, mul_by: Self) -> Self;
    
    fn itanh_mul(self, mul_by: Self) -> Self;
    fn isinh_mul(self, mul_by: Self) -> Self;
    fn icosh_mul(self, mul_by: Self) -> Self;
}

pub trait UintUtils{
    fn usqrt(self) -> Self;
    fn ucbrt(self) -> Self;
    fn usin_mul(self, mul_by: Self) -> Self;
    fn ucos_mul(self, mul_by: Self) -> Self;
    fn utan_mul(self, mul_by: Self) -> Self;
    fn uasin_mul(self, mul_by: Self) -> Self;
    fn uacos_mul(self, mul_by: Self) -> Self;
    fn uatan_mul(self, mul_by: Self) -> Self;
    fn uatan2_mul(self, other: Self, mul_by: Self) -> Self;
    
    fn utanh_mul(self, mul_by: Self) -> Self;
    fn usinh_mul(self, mul_by: Self) -> Self;
    fn ucosh_mul(self, mul_by: Self) -> Self;
}

impl<T: Uint +
        Copy +
        One + Two + Three +
        F32Fmt<F32Fmt = f32> +  
        PartialOrd + 
        Shr<T, Output = T> + 
        Shl<T, Output = T> + 
        Add<T, Output = T> + 
        Mul<T, Output = T>> UintUtils for T {
    #[inline] fn usqrt(self) -> Self {
        if self < Self::TWO { return self; }
        let s_option = Self::usqrt(self >> Self::TWO) << Self::ONE;
        let l_option = s_option + Self::ONE;
        if l_option * l_option > self { s_option } else { l_option }
    }
    #[inline] fn ucbrt(self) -> Self {
        if self < Self::THREE { return self; }
        let s_option = Self::ucbrt(self >> Self::THREE) << Self::ONE;
        let l_option = s_option + Self::ONE;
        if l_option * l_option > self { s_option } else { l_option }
    }
    // note negatives to 0 \/\/\/
    #[inline] fn usin_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().sin() * mul_by.intoF32Fmt()) }
    #[inline] fn ucos_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().cos() * mul_by.intoF32Fmt()) }
    #[inline] fn utan_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().tan() * mul_by.intoF32Fmt()) }
    #[inline] fn uasin_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().asin() * mul_by.intoF32Fmt()) }
    #[inline] fn uacos_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().acos() * mul_by.intoF32Fmt()) }
    #[inline] fn uatan_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().atan() * mul_by.intoF32Fmt()) }
    #[inline] fn uatan2_mul(self, other: Self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().atan2(other.intoF32Fmt()) * mul_by.intoF32Fmt()) }
    
    #[inline] fn usinh_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().sinh() * mul_by.intoF32Fmt()) }
    #[inline] fn ucosh_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().cosh() * mul_by.intoF32Fmt()) }
    #[inline] fn utanh_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().tanh() * mul_by.intoF32Fmt()) }
}

impl<T: Int +
        Copy +
        Zero + One + Two + Three +
        F32Fmt<F32Fmt = f32> + 
        PartialOrd + 
        Shr<T, Output = T> + 
        Shl<T, Output = T> + 
        Add<T, Output = T> + 
        Mul<T, Output = T>> IntUtils for T {
    #[inline] fn isqrt(self) -> Self {
        assert!(self >= Self::ZERO, "sqrt works for only non-negative inputs");
        if self < Self::TWO { return self; }
        let s_option = Self::isqrt(self >> Self::TWO) << Self::ONE;
        let l_option = s_option + Self::ONE;
        if l_option * l_option > self { s_option } else { l_option }
    }
    #[inline] fn icbrt(self) -> Self {
        assert!(self >= Self::ZERO, "sqrt works for only non-negative inputs");
        if self < Self::THREE { return self; }
        let s_option = Self::icbrt(self >> Self::THREE) << Self::ONE;
        let l_option = s_option + Self::ONE;
        if l_option * l_option > self { s_option } else { l_option }
    }
    // note negatives to 0 \/\/\/
    #[inline] fn isin_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().sin() * mul_by.intoF32Fmt()) }
    #[inline] fn icos_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().cos() * mul_by.intoF32Fmt()) }
    #[inline] fn itan_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().tan() * mul_by.intoF32Fmt()) }
    #[inline] fn iasin_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().asin() * mul_by.intoF32Fmt()) }
    #[inline] fn iacos_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().acos() * mul_by.intoF32Fmt()) }
    #[inline] fn iatan_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().atan() * mul_by.intoF32Fmt()) }
    #[inline] fn iatan2_mul(self, other: Self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().atan2(other.intoF32Fmt()) * mul_by.intoF32Fmt()) }
    
    #[inline] fn isinh_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().sinh() * mul_by.intoF32Fmt()) }
    #[inline] fn icosh_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().cosh() * mul_by.intoF32Fmt()) }
    #[inline] fn itanh_mul(self, mul_by: Self) -> Self { Self::fromF32Fmt(self.intoF32Fmt().tanh() * mul_by.intoF32Fmt()) }
}

pub trait Zero { const ZERO: Self; }
impl Zero for bool { const ZERO: Self = false; }
impl Zero for f32 { const ZERO: Self = 0_f32; }
impl Zero for f64 { const ZERO: Self = 0_f64; }
impl Zero for i8    { const ZERO: Self = 0_i8;    } impl Zero for u8    { const ZERO: Self = 0_u8;    }
impl Zero for i16   { const ZERO: Self = 0_i16;   } impl Zero for u16   { const ZERO: Self = 0_u16;   }
impl Zero for i32   { const ZERO: Self = 0_i32;   } impl Zero for u32   { const ZERO: Self = 0_u32;   }
impl Zero for i64   { const ZERO: Self = 0_i64;   } impl Zero for u64   { const ZERO: Self = 0_u64;   }
impl Zero for i128  { const ZERO: Self = 0_i128;  } impl Zero for u128  { const ZERO: Self = 0_u128;  }
impl Zero for isize { const ZERO: Self = 0_isize; } impl Zero for usize { const ZERO: Self = 0_usize; }

pub trait One { const ONE: Self; }
impl One for bool { const ONE: Self = true; }
impl One for f32 { const ONE: Self = 1_f32; }
impl One for f64 { const ONE: Self = 1_f64; }
impl One for i8    { const ONE: Self = 1_i8;    } impl One for u8    { const ONE: Self = 1_u8;    }
impl One for i16   { const ONE: Self = 1_i16;   } impl One for u16   { const ONE: Self = 1_u16;   }
impl One for i32   { const ONE: Self = 1_i32;   } impl One for u32   { const ONE: Self = 1_u32;   }
impl One for i64   { const ONE: Self = 1_i64;   } impl One for u64   { const ONE: Self = 1_u64;   }
impl One for i128  { const ONE: Self = 1_i128;  } impl One for u128  { const ONE: Self = 1_u128;  }
impl One for isize { const ONE: Self = 1_isize; } impl One for usize { const ONE: Self = 1_usize; }

pub trait Two { const TWO: Self; }
impl Two for f32 { const TWO: Self = 2_f32; }
impl Two for f64 { const TWO: Self = 2_f64; }
impl Two for i8    { const TWO: Self = 2_i8;    } impl Two for u8    { const TWO: Self = 2_u8;    }
impl Two for i16   { const TWO: Self = 2_i16;   } impl Two for u16   { const TWO: Self = 2_u16;   }
impl Two for i32   { const TWO: Self = 2_i32;   } impl Two for u32   { const TWO: Self = 2_u32;   }
impl Two for i64   { const TWO: Self = 2_i64;   } impl Two for u64   { const TWO: Self = 2_u64;   }
impl Two for i128  { const TWO: Self = 2_i128;  } impl Two for u128  { const TWO: Self = 2_u128;  }
impl Two for isize { const TWO: Self = 2_isize; } impl Two for usize { const TWO: Self = 2_usize; }

pub trait Three { const THREE: Self; }
impl Three for f32 { const THREE: Self = 3_f32; }
impl Three for f64 { const THREE: Self = 3_f64; }
impl Three for i8    { const THREE: Self = 3_i8;    } impl Three for u8    { const THREE: Self = 3_u8;    }
impl Three for i16   { const THREE: Self = 3_i16;   } impl Three for u16   { const THREE: Self = 3_u16;   }
impl Three for i32   { const THREE: Self = 3_i32;   } impl Three for u32   { const THREE: Self = 3_u32;   }
impl Three for i64   { const THREE: Self = 3_i64;   } impl Three for u64   { const THREE: Self = 3_u64;   }
impl Three for i128  { const THREE: Self = 3_i128;  } impl Three for u128  { const THREE: Self = 3_u128;  }
impl Three for isize { const THREE: Self = 3_isize; } impl Three for usize { const THREE: Self = 3_usize; }

pub trait SignOps { fn ptcopysign(self, sign: Self) -> Self; fn ptsignum(self) -> i8; fn abs(self) -> Self; }
impl SignOps for f32 { #[inline] fn ptcopysign(self, sign: Self) -> Self { self.copysign(sign) } #[inline] fn ptsignum(self) -> i8 { if self == 0.0 { 0_i8 } else if self.signum() > 0.0 { 1_i8 } else { -1_i8 } } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for f64 { #[inline] fn ptcopysign(self, sign: Self) -> Self { self.copysign(sign) } #[inline] fn ptsignum(self) -> i8 { if self == 0.0 { 0_i8 } else if self.signum() > 0.0 { 1_i8 } else { -1_i8 } } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for i8 { #[inline] fn ptcopysign(self, sign: Self) -> Self { if sign >> 7 ^ self >> 7 != 0 { -self } else { self } } #[inline] fn ptsignum(self) -> i8 { self.signum() } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for i16 { #[inline] fn ptcopysign(self, sign: Self) -> Self { if sign >> 15 ^ self >> 15 != 0  { -self } else { self } } #[inline] fn ptsignum(self) -> i8 { self.signum() as i8 } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for i32 { #[inline] fn ptcopysign(self, sign: Self) -> Self { if sign >> 31 ^ self >> 31 != 0  { -self } else { self } } #[inline] fn ptsignum(self) -> i8 { self.signum() as i8 } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for i64 { #[inline] fn ptcopysign(self, sign: Self) -> Self { if sign >> 63 ^ self >> 63 != 0  { -self } else { self } } #[inline] fn ptsignum(self) -> i8 { self.signum() as i8 } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for i128 { #[inline] fn ptcopysign(self, sign: Self) -> Self { if sign >> 127 ^ self >> 127 != 0  { -self } else { self } } #[inline] fn ptsignum(self) -> i8 { self.signum() as i8 } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for isize { #[inline] fn ptcopysign(self, sign: Self) -> Self { if sign.reverse_bits() >> 1 ^ self.reverse_bits() >> 1 != 0  { -self } else { self } } #[inline] fn ptsignum(self) -> i8 { self.signum() as i8 } #[inline] fn abs(self) -> Self { self.abs() } }
impl SignOps for u8 { #[inline] fn ptcopysign(self, _: Self) -> Self { self } #[inline] fn ptsignum(self) -> i8 { if self != 0_u8 { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { self } }
impl SignOps for u16 { #[inline] fn ptcopysign(self, _: Self) -> Self { self } #[inline] fn ptsignum(self) -> i8 { if self != 0_u16 { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { self } }
impl SignOps for u32 { #[inline] fn ptcopysign(self, _: Self) -> Self { self } #[inline] fn ptsignum(self) -> i8 { if self != 0_u32 { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { self } }
impl SignOps for u64 { #[inline] fn ptcopysign(self, _: Self) -> Self { self } #[inline] fn ptsignum(self) -> i8 { if self != 0_u64 { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { self } }
impl SignOps for u128 { #[inline] fn ptcopysign(self, _: Self) -> Self { self } #[inline] fn ptsignum(self) -> i8 { if self != 0_u128 { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { self } }
impl SignOps for usize { #[inline] fn ptcopysign(self, _: Self) -> Self { self } #[inline] fn ptsignum(self) -> i8 { if self != 0_usize { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { self } }
impl SignOps for bool { #[inline] fn ptcopysign(self, sign: Self) -> Self { sign } #[inline] fn ptsignum(self) -> i8 { if self { 1_i8 } else { 0_i8 } } #[inline] fn abs(self) -> Self { false } }

pub trait Primitive {}
impl Primitive for f32 {}
impl Primitive for f64 {}
impl Primitive for i8 {}
impl Primitive for i16 {}
impl Primitive for i32 {}
impl Primitive for i64 {}
impl Primitive for i128 {}
impl Primitive for isize {}
impl Primitive for u8 {}
impl Primitive for u16 {}
impl Primitive for u32 {}
impl Primitive for u64 {}
impl Primitive for u128 {}
impl Primitive for usize {}
impl Primitive for bool {}

pub trait Float {}
impl Float for f32 {}
impl Float for f64 {}

pub trait Int {}
impl Int for i8 {}
impl Int for i16 {}
impl Int for i32 {}
impl Int for i64 {}
impl Int for i128 {}
impl Int for isize {}

pub trait Uint {}
impl Uint for u8 {}
impl Uint for u16 {}
impl Uint for u32 {}
impl Uint for u64 {}
impl Uint for u128 {}
impl Uint for usize {}
impl Uint for bool {}

pub trait F32Fmt {
    type F32Fmt: F32Fmt + Copy + Clone + PartialEq
        + Add<Self::F32Fmt, Output = Self::F32Fmt>
        + Sub<Self::F32Fmt, Output = Self::F32Fmt>
        + Mul<Self::F32Fmt, Output = Self::F32Fmt>
        + Div<Self::F32Fmt, Output = Self::F32Fmt>
        + Rem<Self::F32Fmt, Output = Self::F32Fmt>
        + Neg<Output = Self::F32Fmt> + SignOps
        + Zero + One + Two;
    #[allow(non_snake_case)]
    fn intoF32Fmt(self) -> Self::F32Fmt;
    #[allow(non_snake_case)]
    fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self;
    fn sqrt(self) -> Self;
    fn cbrt(self) -> Self;
    fn f32_const_mul(self, constant: f32) -> Self;
    fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    
    fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
    fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized;
}
impl F32Fmt for f32 { type F32Fmt = Self;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self } 
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt }
    #[inline] fn sqrt(self) -> Self { Self::sqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::cbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { self * constant }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::sin(self) * mul_by }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::cos(self) * mul_by }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::tan(self) * mul_by }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::asin(self) * mul_by }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::acos(self) * mul_by }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::atan(self) * mul_by }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::atan2(self, other) * mul_by }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::sinh(self) * mul_by }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::cosh(self) * mul_by }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::tanh(self) * mul_by }
}
impl F32Fmt for f64 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::sqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::cbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self where Self: Mul<f64, Output = Self> { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::sin(self) * mul_by }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::cos(self) * mul_by }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::tan(self) * mul_by }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::asin(self) * mul_by }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::acos(self) * mul_by }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::atan(self) * mul_by }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::atan2(self, other) * mul_by }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::sinh(self) * mul_by }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::cosh(self) * mul_by }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::tanh(self) * mul_by }
}
impl F32Fmt for i8 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::isqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::icbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itanh_mul(self, mul_by) }
}
impl F32Fmt for i16 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 } 
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::isqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::icbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itanh_mul(self, mul_by) }
}
impl F32Fmt for i32 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::isqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::icbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itanh_mul(self, mul_by) }
}
impl F32Fmt for i64 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::isqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::icbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itanh_mul(self, mul_by) }
}
impl F32Fmt for i128 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::isqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::icbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itanh_mul(self, mul_by) }
}
impl F32Fmt for isize {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::isqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::icbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::iatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::isinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::icosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::itanh_mul(self, mul_by) }
}
impl F32Fmt for u8 { 
    type F32Fmt = f32; 
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 } 
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::usqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::ucbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utanh_mul(self, mul_by) }
}
impl F32Fmt for u16 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::usqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::ucbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utanh_mul(self, mul_by) }
}
impl F32Fmt for u32 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::usqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::ucbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utanh_mul(self, mul_by) }
}
impl F32Fmt for u64 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::usqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::ucbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utanh_mul(self, mul_by) }
}
impl F32Fmt for u128 {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::usqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::ucbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utanh_mul(self, mul_by) }
}
impl F32Fmt for usize {
    type F32Fmt = f32;
    #[inline] fn intoF32Fmt(self) -> Self::F32Fmt { self as f32 }
    #[inline] fn fromF32Fmt(f32_fmt: Self::F32Fmt) -> Self { f32_fmt as Self }
    #[inline] fn sqrt(self) -> Self { Self::usqrt(self) }
    #[inline] fn cbrt(self) -> Self { Self::ucbrt(self) }
    #[inline] fn f32_const_mul(self, constant: f32) -> Self { (self as f32 * constant) as Self }
    #[inline] fn sin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usin_mul(self, mul_by) }
    #[inline] fn cos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucos_mul(self, mul_by) }
    #[inline] fn tan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utan_mul(self, mul_by) }
    #[inline] fn asin_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uasin_mul(self, mul_by) }
    #[inline] fn acos_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uacos_mul(self, mul_by) }
    #[inline] fn atan_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan_mul(self, mul_by) }
    #[inline] fn atan2_mul(self, other: Self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::uatan2_mul(self, other, mul_by) }
    
    #[inline] fn sinh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::usinh_mul(self, mul_by) }
    #[inline] fn cosh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::ucosh_mul(self, mul_by) }
    #[inline] fn tanh_mul(self, mul_by: Self) -> Self where Self: Mul<Self, Output = Self> + Sized { Self::utanh_mul(self, mul_by) }
}
