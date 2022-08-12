use std::ops::Mul;

use macroquad::prelude::{Vec2, };

pub struct Chunk {
    w: usize,
    points: Vec<usize>,
    center: Vec2,
    area: usize,
    // The major axis length
    major: f32,
    // The minor axis length
    minor: f32,
    // The direction of the major axis.
    axis: Vec2,
}

impl Chunk {
    pub fn find(w: usize, mut pixels: Vec<bool>) -> Vec<Self> {
        let mut out = Vec::new();
        while let Some(points) = contiguous_pixels(w, &mut pixels) {
            out.push(Chunk::new(w, points));
        }
        out
    }
    pub fn new(w: usize, points: Vec<usize>) -> Self {
        let area = points.len();
        let mut center = Vec2::ZERO;
        for p in points.iter().copied() {
            let x = (p % w) as f32;
            let y = (p / w) as f32;
            center += Vec2::new(x, y);
        }
        center /= area as f32;
        let mut x2 = 0.0;
        let mut y2 = 0.0;
        let mut xy = 0.0;
        for p in points.iter().copied() {
            let dx = (p % w) as f32 - center.x;
            let dy = (p / w) as f32 - center.y;
            x2 += dx * dx;
            y2 += dy * dy;
            xy += dx * dy;
        }
        // v = (a b)
        // x2*a + xy*b = e*a
        // xy*a + y2*b = e*b
        //
        // x2*ab + xy*b**2 = e*ab
        // xy*a**2 + y2*ab = e*ab
        //
        // xy*a**2 + (x2-y2)ab - xy*b**2 = 0
        //
        // xy*(a/b)**2 + (x2-y2)*(a/b) - xy = 0
        //
        // c == (a/b)
        //
        // c = ((y2-x2) +/- sqrt((x2-y2)**2 + 4*xy**2)) / (2*xy)
        let ax1 = Vec2::new(
            ((y2 - x2) + ((x2 - y2).powi(2) + 4.0 * xy.powi(2)).sqrt()) / (2.0 * xy),
            1.0,
        )
        .normalize();
        let ax2 = Vec2::new(
            ((y2 - x2) - ((x2 - y2).powi(2) + 4.0 * xy.powi(2)).sqrt()) / (2.0 * xy),
            1.0,
        )
        .normalize();
        // e = x2 + xy*(b/a) = y2 + xy*(a/b)
        let e1 = if ax1.x.abs() > ax1.y.abs() {
            y2 + xy * ax1.x / ax1.y
        } else {
            x2 + xy * ax1.y / ax1.x
        };
        let e2 = if ax2.x.abs() > ax2.y.abs() {
            y2 + xy * ax2.x / ax2.y
        } else {
            x2 + xy * ax2.y / ax2.x
        };
        let (major, minor, axis) = if xy == 0.0 {
            if x2 > y2 {
                (x2.sqrt(), y2.sqrt(), Vec2::new(1.0, 0.0))
            } else {
                (y2.sqrt(), x2.sqrt(), Vec2::new(0.0, 1.0))
            }
        } else if e1 > e2 {
            (e1.sqrt(), e2.sqrt(), ax1)
        } else {
            (e2.sqrt(), e1.sqrt(), ax2)
        };
        Chunk {
            w,
            area,
            points,
            center,
            major,
            minor,
            axis,
        }
    }
}

fn contiguous_pixels(w: usize, pixels: &mut [bool]) -> Option<Vec<usize>> {
    let mut out = Vec::new();

    let (p, _) = pixels.iter().enumerate().filter(|(_, p)| **p).next()?;
    pixels[p] = false;
    let mut todo = vec![p];
    while let Some(p) = todo.pop() {
        out.push(p);
        let x = p % w;
        let y = p / w;
        if x > 0 && pixels[p - 1] {
            todo.push(p - 1);
            pixels[p - 1] = false;
        }
        if x < w - 1 && pixels[p + 1] {
            todo.push(p + 1);
            pixels[p + 1] = false;
        }
        if y > 0 && pixels[p - w] {
            todo.push(p - w);
            pixels[p - w] = false;
        }
        if y < w - 1 && pixels[p + w] {
            todo.push(p + w);
            pixels[p + w] = false;
        }
    }
    Some(out)
}

#[test]
fn chunks_test() {
    let pixels = vec![false, false, true, false, false, true, false, false, true];
    let w = 3;
    let chunks = Chunk::find(w, pixels);
    assert_eq!(1, chunks.len());
    assert_eq!(3, chunks[0].area);
    assert_eq!(Vec2::new(2.0, 1.0), chunks[0].center);
    assert_eq!(Vec2::new(0.0, 1.0), chunks[0].axis);
    assert_eq!(2.0_f32.sqrt(), chunks[0].major);
    assert_eq!(0.0, chunks[0].minor);

    let pixels = vec![false, false, true, true, false, true, true, false, true];
    let w = 3;
    let chunks = Chunk::find(w, pixels);
    assert_eq!(2, chunks.len());
    assert_eq!(3, chunks[0].area);
    assert_eq!(2, chunks[1].area);
    assert_eq!(Vec2::new(2.0, 1.0), chunks[0].center);
    assert_eq!(Vec2::new(0.0, 1.0), chunks[0].axis);
    assert_eq!(2.0_f32.sqrt(), chunks[0].major);
    assert_eq!(0.0, chunks[0].minor);
}


#[derive(Copy, Clone)]
pub struct Transform {
    full_angle: f32,
    sin: f32,
    cos: f32,
    scale_major: f32,
    major_axis: Vec2,
    scale_minor: f32,
    minor_axis: Vec2,
    center: Vec2,
    translation: Vec2,
}

impl Transform {
    pub fn new(o: &Chunk, n: &Chunk) -> Self {
        let full_angle = if o.axis.dot(n.axis) > 0.0 {
            n.axis.angle_between(o.axis)
        } else {
            n.axis.angle_between(-o.axis)
        };
        let (sin, cos) = full_angle.sin_cos();
        let major_axis = o.axis;
        let minor_axis = o.axis.perp();
        let scale_major = if o.major > 0.0 && n.major > 0.0 {
            n.major / o.major
        } else {
            1.0
        };
        let scale_minor = if o.minor > 0.0 && n.minor > 0.0 {
            n.minor / o.minor
        } else {
            1.0
        };
        let center = o.center;
        let translation = n.center - o.center;
        Transform { full_angle, sin, cos, scale_major, major_axis, scale_minor, minor_axis, center, translation }
    }
    pub fn reverse(mut self) -> Self {
        self.full_angle *= -1.0;
        self.major_axis = Vec2::new(self.major_axis.x*self.cos + self.major_axis.y*self.sin, self.major_axis.y*self.cos - self.major_axis.x*self.sin);
        self.minor_axis = Vec2::new(self.minor_axis.x*self.cos + self.minor_axis.y*self.sin, self.minor_axis.y*self.cos - self.minor_axis.x*self.sin);
        self.sin *= -1.0;
        self.center = self.center + self.translation;
        self.translation *= -1.0;
        self.scale_major = 1.0/self.scale_major;
        self.scale_minor = 1.0/self.scale_minor;

        self
    }
    pub fn scale(mut self, f: f32) -> Self {
        self.full_angle *= f;
        (self.sin, self.cos) = self.full_angle.sin_cos();
        self.scale_major = self.scale_major.powf(f);
        self.scale_minor = self.scale_minor.powf(f);
        self.translation *= f;
        self
    }
}

impl Mul<Transform> for f32 {
    type Output = Transform;
    fn mul(self, rhs: Transform) -> Self::Output {
        rhs.scale(self)
    }
}

impl Mul<Vec2> for Transform {
    type Output = Vec2;

    fn mul(self, rhs: Vec2) -> Self::Output {
        let rhs = rhs - self.center;
        let v = self.major_axis*self.major_axis.dot(rhs)*self.scale_major + self.minor_axis*self.minor_axis.dot(rhs)*self.scale_minor;
        self.center + self.translation + Vec2::new(v.x*self.cos + v.y*self.sin, v.y*self.cos - v.x*self.sin)
    }
}