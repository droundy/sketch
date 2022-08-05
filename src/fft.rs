// use macroquad::texture::Image;

// use num_complex::Complex;

// pub fn fft(img: &Image) -> Vec<Complex<f64>> {
//     let mut out: Vec<_> = img
//         .get_image_data()
//         .iter()
//         .map(|p| Complex::new(p[3] as f64 / 255.0, 0.0))
//         .collect();
//     assert_eq!(out.len(), img.width as usize * img.height as usize);
//     fft2d::slice::fft_2d(img.width as usize, img.height as usize, &mut out);
//     out
// }

// pub fn ifft(img: &mut Image, mut data: Vec<Complex<f64>>, num_pixels: usize) {
//     fft2d::slice::ifft_2d(img.width as usize, img.height as usize, &mut data);
//     let mut lo = 0.0;
//     let mut hi = data.len() as f64 * 100.0;
//     while hi - lo > 0.001 {
//         let mid = 0.5 * (hi + lo);
//         let count = data.iter().filter(|v| v.re > mid).count();
//         if count > num_pixels {
//             lo = mid;
//         } else {
//             hi = mid;
//         }
//     }

//     for (i, v) in img.get_image_data_mut().iter_mut().enumerate() {
//         if data[i].re > lo {
//             *v = [255; 4];
//         } else {
//             *v = [0; 4];
//         }
//     }
// }

// #[test]
// fn fft_inverse_itself() {
//     let width = 1024;
//     let height = 768;
//     let mut v = vec![Complex::new(0.0, 0.0); width * height];
//     v[1].re = 5.0;
//     fft2d::slice::fft_2d(width, height, &mut v);
//     fft2d::slice::ifft_2d(width, height, &mut v);
//     assert_eq!(v[0].re, 0.0);
//     assert_eq!(v[2].re, 0.0);
//     assert!(v[1].re > 0.0);
//     for i in 2.. width * height {
//         println!(" {i}   {:?}", v[i]);
//         assert!(v[i].re.abs() < 1e-7);
//     }
// }

use std::f32::consts::PI;
use std::sync::Arc;

use rustfft::num_complex::Complex;
use rustfft::{Fft, FftPlanner};
// use rustfft::num_traits::Zero;

use crate::layer::clamp;

pub struct Ffts {
    height: usize,
    width: usize,

    for_width: Arc<dyn Fft<f32>>,
    bac_width: Arc<dyn Fft<f32>>,
    for_height: Arc<dyn Fft<f32>>,
    bac_height: Arc<dyn Fft<f32>>,
}

impl Ffts {
    pub fn fft(&mut self, img: &[[u8; 4]]) -> Vec<Complex<f32>> {
        assert_eq!(img.len(), self.width * self.height);
        let mut output: Vec<_> = vec![Complex::new(0.0, 0.0); self.width * self.height];
        for i in 0..self.height {
            for j in 0..self.width {
                let idx = i * self.width + j;
                output[idx] = Complex::new(img[idx][3] as f32 / 255.0, 0.0);
            }
            self.for_width
                .process(&mut output[i * self.width..(i + 1) * self.width]);
        }
        let mut col = vec![Complex::new(0.0, 0.0); self.height];
        for j in 0..self.width {
            for i in 0..self.height {
                col[i] = output[i * self.width + j];
            }
            self.for_height.process(&mut col);
            for i in 0..self.height {
                output[i * self.width + j] = col[i];
            }
        }
        // Do something evil and store ln(fft) here instead.
        for i in 0..self.height {
            for j in 0..self.width {
                let idx = i * self.width + j;
                output[idx] = output[idx].ln();
            }
        }
        for i in 1..self.height {
            while output[i * self.width].im < output[(i - 1) * self.width].im - PI {
                output[i * self.width].im += 2.0 * PI;
            }
            while output[i * self.width].im > output[(i - 1) * self.width].im + PI {
                output[i * self.width].im -= 2.0 * PI;
            }
        }
        for i in 0..self.height {
            for j in 1..self.width {
                let idx = i * self.width + j;
                while output[idx].im < output[idx - 1].im - PI {
                    output[idx].im += 2.0 * PI;
                }
                while output[idx].im > output[idx - 1].im + PI {
                    output[idx].im -= 2.0 * PI;
                }
            }
        }
        let mut peak = -1e20;
        let mut peak_phase = 0.0;
        for v in output.iter() {
            if v.re > peak {
                peak = v.re;
                peak_phase = v.im;
            }
        }
        for v in output.iter_mut() {
            v.im -= peak_phase;
        }
        output
    }
    pub fn ifft(&mut self, img: &mut [[u8; 4]], mut data: Vec<Complex<f32>>) {
        assert_eq!(img.len(), self.width * self.height);
        assert_eq!(data.len(), img.len());
        let mut col = vec![Complex::new(0.0, 0.0); self.height];
        for j in 0..self.width {
            for i in 0..self.height {
                col[i] = data[i * self.width + j];
            }
            self.bac_height.process(&mut col);
            for i in 0..self.height {
                data[i * self.width + j] = col[i];
            }
        }
        for i in 0..self.height {
            self.bac_width
                .process(&mut data[i * self.width..(i + 1) * self.width]);
            for j in 0..self.width {
                let idx = i * self.width + j;
                img[idx] = [
                    255,
                    255,
                    255,
                    clamp(
                        0.0,
                        255.0,
                        255.0 * data[idx].re / self.width as f32 / self.height as f32,
                    )
                    .round() as u8,
                ];
            }
        }
    }

    pub fn new(width: usize, height: usize) -> Self {
        // make a planner
        let mut real_planner = FftPlanner::<f32>::new();
        Ffts {
            height,
            width,
            for_width: real_planner.plan_fft_forward(width),
            bac_width: real_planner.plan_fft_inverse(width),
            for_height: real_planner.plan_fft_forward(height),
            bac_height: real_planner.plan_fft_inverse(height),
        }
    }
}
