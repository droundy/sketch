use std::ops::Mul;

use glam::Vec2;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use tinyset::SetUsize;

/// A connection between two keyframes.
#[derive(Clone, Serialize, Deserialize)]
pub struct Tween {
    chunks: Vec<ChunkTween>,
}

impl Tween {
    pub fn new(w: usize, before: SetUsize, after: SetUsize) -> Self {
        let mut chunks = Vec::new();

        let mut before_chunks = Chunk::find(w, before);
        before_chunks.sort_by(|a, b| b.area.cmp(&a.area));
        let mut after_chunks = Chunk::find(w, after);
        after_chunks.sort_by(|a, b| b.area.cmp(&a.area));

        let nchunks = std::cmp::min(before_chunks.len(), after_chunks.len());
        if nchunks == 0 {
            return Tween { chunks: Vec::new() };
        }
        // FIXME this is a hokwy way to pair up the chunks.
        for _ in 0..1000 {
            let i = rand::random::<usize>() % nchunks;
            let j = rand::random::<usize>() % nchunks;
            let v_before = before_chunks[i].center - before_chunks[j].center;
            let v_after = after_chunks[i].center - after_chunks[j].center;
            if v_after.dot(v_before) < 0.0 {
                after_chunks.swap(i, j);
            }
        }
        for (before, after) in before_chunks.into_iter().zip(after_chunks.into_iter()) {
            chunks.push(ChunkTween::new(w, before, after));
        }
        Tween { chunks }
    }

    /// Draw this Tween to the pixel buffer
    pub fn draw(&self, fraction: f32, pixels: &mut [bool]) {
        for c in self.chunks.iter() {
            c.draw(fraction, pixels);
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct ChunkTween {
    w: usize,
    connections: Vec<(usize, usize)>,
    transform: Transform,
}

impl ChunkTween {
    /// Figure out how to tween between a before and after [`Chunk`].
    ///
    /// FIXME: Want to respect transform in connecting pixels!
    fn new(w: usize, before: Chunk, after: Chunk) -> Self {
        let transform = Transform::new(&before, &after);

        let mut connections = Vec::new();

        let mut before_pixels = before.points.clone();
        let before_outline = outline(w, &mut before_pixels).unwrap();
        let before_strokes = strokes_from_outline(w, &before_outline);
        println!("strokes before are {}", before_strokes.len());
        let mut after_pixels = after.points.clone();
        let after_outline = outline(w, &mut after_pixels).unwrap();
        let after_strokes = strokes_from_outline(w, &after_outline);
        println!("strokes after are {}", after_strokes.len());
        // for s in after_strokes.iter() {
        //     println!("    {s:?}")
        // }

        let mut i = 0;
        let mut j = 0;
        while i < before_outline.len() && j < after_outline.len() {
            connections.push((before_outline[i], after_outline[j]));
            let ifrac = (i + 1) as f64 / before_outline.len() as f64;
            let jfrac = (j + 1) as f64 / after_outline.len() as f64;
            if ifrac < jfrac {
                i += 1;
            } else {
                j += 1;
            }
        }
        let mut todo = connections.clone();

        while !before_pixels.is_empty() || !after_pixels.is_empty() {
            let mut more = Vec::new();
            for (b, a) in todo.drain(..) {
                if before_pixels.contains(b + 1) {
                    if after_pixels.contains(a + 1) {
                        more.push((b + 1, a + 1));
                    } else {
                        more.push((b + 1, a));
                    }
                } else if after_pixels.contains(a + 1) {
                    more.push((b, a + 1));
                }
                if before_pixels.contains(b + w) {
                    if after_pixels.contains(a + w) {
                        more.push((b + w, a + w));
                    } else {
                        more.push((b + w, a));
                    }
                } else if after_pixels.contains(a + w) {
                    more.push((b, a + w));
                }
                if before_pixels.contains(b.wrapping_sub(1)) {
                    if after_pixels.contains(a.wrapping_sub(1)) {
                        more.push((b.wrapping_sub(1), a.wrapping_sub(1)));
                    } else {
                        more.push((b.wrapping_sub(1), a));
                    }
                } else if after_pixels.contains(a.wrapping_sub(1)) {
                    more.push((b, a.wrapping_sub(1)));
                }
                if before_pixels.contains(b.wrapping_sub(w)) {
                    if after_pixels.contains(a.wrapping_sub(w)) {
                        more.push((b.wrapping_sub(w), a.wrapping_sub(w)));
                    } else {
                        more.push((b.wrapping_sub(w), a));
                    }
                } else if after_pixels.contains(a.wrapping_sub(w)) {
                    more.push((b, a.wrapping_sub(w)));
                }
            }
            more.sort();
            more.dedup();
            for (b, a) in more.iter().copied() {
                before_pixels.remove(b);
                after_pixels.remove(a);
            }
            connections.extend(more.iter().copied());
            todo = more;
        }

        connections.sort();
        connections.dedup();

        ChunkTween {
            w,
            connections,
            transform,
        }
    }

    /// Draw this ChunkTween to the pixel buffer
    fn draw(&self, fraction: f32, pixels: &mut [bool]) {
        assert!(fraction >= 0.0);
        assert!(fraction <= 1.0);
        let reverse_transform = (1.0 - fraction) * self.transform.reverse();
        let transform = fraction * self.transform;
        for &(b, a) in self.connections.iter() {
            let b = transform * Vec2::new((b % self.w) as f32, (b / self.w) as f32);
            let a = reverse_transform * Vec2::new((a % self.w) as f32, (a / self.w) as f32);
            let p = fraction * a + (1.0 - fraction) * b;
            let w = self.w;
            let idx0 = p.x as usize + (p.y as usize) * w;
            for idx in [idx0, idx0 + 1, idx0 + w, idx0 + w + 1] {
                if let Some(p) = pixels.get_mut(idx) {
                    *p = true;
                }
            }
        }
    }
}

pub struct Chunk {
    points: SetUsize,
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
    pub fn find(w: usize, mut pixels: SetUsize) -> Vec<Self> {
        let mut out = Vec::new();
        while let Some(points) = contiguous_pixels(w, &mut pixels) {
            out.push(Chunk::new(w, points));
        }
        out
    }
    pub fn new(w: usize, points: SetUsize) -> Self {
        let area = points.len();
        let mut center = Vec2::ZERO;
        for p in points.iter() {
            let x = (p % w) as f32;
            let y = (p / w) as f32;
            center += Vec2::new(x, y);
        }
        center /= area as f32;
        let mut x2 = 0.0;
        let mut y2 = 0.0;
        let mut xy = 0.0;
        for p in points.iter() {
            let dx = (p % w) as f32 - center.x;
            let dy = (p / w) as f32 - center.y;
            x2 += dx * dx;
            y2 += dy * dy;
            xy += dx * dy;
        }
        x2 /= points.len() as f32;
        y2 /= points.len() as f32;
        xy /= points.len() as f32;
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
        let ax1 = if xy == 0.0 {
            Vec2::new(1.0, 0.0)
        } else {
            Vec2::new(
                (-(y2 - x2) - ((x2 - y2).powi(2) + 4.0 * xy.powi(2)).sqrt()) / (2.0 * xy),
                1.0,
            )
            .normalize()
        };
        let ax2 = if xy == 0.0 {
            Vec2::new(0.0, 1.0)
        } else {
            Vec2::new(
                -((y2 - x2) - ((x2 - y2).powi(2) + 4.0 * xy.powi(2)).sqrt()) / (2.0 * xy),
                1.0,
            )
            .normalize()
        };
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
        // println!("axis 1: {} {}  gives {e1}", ax1.x, ax1.y);
        // println!("axis 2: {} {}  gives {e2}", ax2.x, ax2.y);
        Chunk {
            area,
            points,
            center,
            major,
            minor,
            axis,
        }
    }
}

fn contiguous_pixels(w: usize, pixels: &mut SetUsize) -> Option<SetUsize> {
    let mut out = SetUsize::new();

    let p = pixels.iter().next()?;
    pixels.remove(p);
    let mut todo = vec![p];
    while let Some(p) = todo.pop() {
        out.insert(p);
        if p > 0 && pixels.contains(p - 1) {
            todo.push(p - 1);
            pixels.remove(p - 1);
        }
        if pixels.contains(p + 1) {
            todo.push(p + 1);
            pixels.remove(p + 1);
        }
        if p >= w && pixels.contains(p - w) {
            todo.push(p - w);
            pixels.remove(p - w);
        }
        if pixels.contains(p + w) {
            todo.push(p + w);
            pixels.remove(p + w);
        }
    }
    Some(out)
}

#[derive(Clone, Copy)]
struct Width(usize);
impl Width {
    fn dist_sqr(self, a: usize, b: usize) -> f32 {
        ((a % self.0) as f32 - (b % self.0) as f32).powi(2)
            + ((a / self.0) as f32 - (b / self.0) as f32).powi(2)
    }
}

fn closest_local_minimum(w: Width, outline: &[usize], i: usize) -> Option<usize> {
    let mut min_dist2 = None;
    let mut min_index = 0;
    for j in 0..outline.len() {
        let dj = w.dist_sqr(outline[i], outline[j]);
        let djp1 = w.dist_sqr(outline[i], outline[(j + 1) % outline.len()]);
        let djm1 = w.dist_sqr(outline[i], outline[(j + outline.len() - 1) % outline.len()]);
        if j != i && dj < djp1 && dj < djm1 {
            // This is a local minimum.
            if let Some(existing_min) = min_dist2 {
                if dj < existing_min {
                    min_dist2 = Some(dj);
                    min_index = j;
                }
            } else {
                min_dist2 = Some(dj);
                min_index = j;
            }
        }
    }
    if min_dist2.is_some() {
        Some(min_index)
    } else {
        None
    }
}

fn mutual_minimum(w: Width, outline: &[usize], i: usize) -> Option<usize> {
    let j = closest_local_minimum(w, outline, i)?;
    if closest_local_minimum(w, outline, j)? != i {
        None
    } else {
        Some(j)
    }
}

fn strokes_from_outline(w: usize, outline: &[usize]) -> Vec<Vec<(usize, usize)>> {
    let w = Width(w);
    let mut strokes = Vec::new();
    let mut possibilities = SetUsize::from_iter(0..outline.len());
    // while let Some(p) = possibilities.iter().next() {
    for i in 0..outline.len() {
        if possibilities.contains(i) {
            if let Some(mut j) = mutual_minimum(w, outline, i) {
                // We have found a new stroke, defined as a point that has
                // another point closest to it which is also a local minimum.
                possibilities.remove(j);
                let mut stroke = vec![(i, j)];
                // Now we want to find other members of this stroke.
                let mut left = i;
                let mut right = j;
                let mut previous_lr = (0, 0);
                while previous_lr != (left, right) {
                    previous_lr = (left, right);
                    let l = (left + 1) % outline.len();
                    possibilities.remove(l);
                    if let Some(r) = closest_local_minimum(w, outline, l) {
                        possibilities.remove(r);
                        if r == right {
                            left = l;
                            stroke.push((left, right));
                        } else if r == (right + outline.len() - 1) % outline.len() {
                            left = l;
                            right = r;
                            stroke.push((left, right));
                        } else if r == (right + 2 * outline.len() - 2) % outline.len() {
                            let intermediate = (right + outline.len() - 1) % outline.len();
                            possibilities.remove(intermediate);
                            stroke.push((left, intermediate));
                            left = l;
                            stroke.push((left, intermediate));
                            right = r;
                            stroke.push((left, right));
                        }
                    }
                    let r = (right + outline.len() - 1) % outline.len();
                    possibilities.remove(r);
                    if let Some(l) = closest_local_minimum(w, outline, r) {
                        possibilities.remove(l);
                        if l == left {
                            right = r;
                            stroke.push((left, right));
                        } else if l == (left + 1) % outline.len() {
                            right = r;
                            left = l;
                            stroke.push((left, right));
                        } else if l == (left + 2) % outline.len() {
                            let intermediate = (left + 1) % outline.len();
                            possibilities.remove(intermediate);
                            stroke.push((intermediate, right));
                            right = r;
                            stroke.push((intermediate, right));
                            left = l;
                            stroke.push((left, right));
                        }
                    }
                }
                stroke.sort_unstable();
                strokes.push(stroke);
            }
        }
    }
    strokes.retain(|v| v.len() > 1);
    strokes
    // strokes
    //     .into_iter()
    //     .map(|v| {
    //         v.into_iter()
    //             .map(|(i, j)| (outline[i], outline[j]))
    //             .collect()
    //     })
    //     .collect()
}

fn outline(w: usize, pixels: &mut SetUsize) -> Option<Vec<usize>> {
    let mut out = Vec::new();

    let mut start = pixels.iter().next()?;
    let mut best_diag = (start % w) + (start / w);
    for p in pixels.iter() {
        if (p % w) + (p / w) < best_diag {
            start = p;
            best_diag = (p % w) + (p / w);
        }
    }
    out.push(start);
    let mut last = start;
    let mut next = if pixels.contains(last + 1) {
        last + 1
    } else if pixels.contains(last + w) {
        last + w
    } else {
        return Some(out);
    };
    pixels.remove(next);
    out.push(next);
    loop {
        let n = if next == last + 1 {
            if pixels.contains(next.wrapping_sub(w)) {
                next.wrapping_sub(w)
            } else if pixels.contains(next + 1) {
                next + 1
            } else if pixels.contains(next + w) {
                next + w
            } else {
                break;
            }
        } else if next == last + w {
            if pixels.contains(next + 1) {
                next + 1
            } else if pixels.contains(next + w) {
                next + w
            } else if pixels.contains(next.wrapping_sub(1)) {
                next.wrapping_sub(1)
            } else {
                break;
            }
        } else if next == last.wrapping_sub(1) {
            if pixels.contains(next + w) {
                next + w
            } else if pixels.contains(next.wrapping_sub(1)) {
                next.wrapping_sub(1)
            } else if pixels.contains(next.wrapping_sub(w)) {
                next.wrapping_sub(w)
            } else {
                break;
            }
        } else if next == last.wrapping_sub(w) {
            if pixels.contains(next.wrapping_sub(1)) {
                next.wrapping_sub(1)
            } else if pixels.contains(next.wrapping_sub(w)) {
                next.wrapping_sub(w)
            } else if pixels.contains(next + 1) {
                next + 1
            } else {
                break;
            }
        } else {
            unreachable!()
        };
        pixels.remove(n);
        out.push(n);
        last = next;
        next = n;
        if n == start {
            break;
        }
    }
    Some(out)
}

#[test]
fn chunks_test() {
    let pixels = SetUsize::from_iter([2, 5, 8]);
    let w = 3;
    let chunks = Chunk::find(w, pixels);
    assert_eq!(1, chunks.len());
    assert_eq!(3, chunks[0].area);
    assert_eq!(Vec2::new(2.0, 1.0), chunks[0].center);
    assert_eq!(Vec2::new(0.0, 1.0), chunks[0].axis);
    assert_eq!(2.0_f32.sqrt(), chunks[0].major);
    assert_eq!(0.0, chunks[0].minor);

    let pixels = SetUsize::from_iter([2, 3, 5, 6, 8]);
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

#[derive(Copy, Clone, Serialize, Deserialize)]
pub struct Transform {
    full_angle: f32,
    sin: f32,
    cos: f32,
    scale_major: f32,
    major_axis: Vec2,
    scale_minor: f32,
    center: Vec2,
    translation: Vec2,
}

impl Transform {
    pub fn new(o: &Chunk, n: &Chunk) -> Self {
        let mut full_angle = if o.axis.dot(n.axis) > 0.0 {
            -n.axis.angle_between(o.axis)
        } else {
            -n.axis.angle_between(-o.axis)
        };
        full_angle = 0.0;
        // println!("full_angle is {full_angle}");
        let (sin, cos) = full_angle.sin_cos();
        let major_axis = o.axis;
        let scale_major = if false && o.major > 0.0 && n.major > 0.0 {
            n.major / o.major
        } else {
            1.0
        };
        let scale_minor = if false && o.minor > 0.0 && n.minor > 0.0 {
            n.minor / o.minor
        } else {
            1.0
        };
        if o.major / o.minor < 1.1 || n.major / n.minor < 1.1 {
            full_angle = 0.0;
        }
        // println!("scales are: {scale_major} and {scale_minor}");
        // println!("old major {} and minor {}", o.major, o.minor);
        // println!("new major {} and minor {}", n.major, n.minor);
        // println!("the major axis is {} {}", major_axis.x, major_axis.y);
        let center = o.center;
        let translation = n.center - o.center;
        Transform {
            full_angle,
            sin,
            cos,
            scale_major,
            major_axis,
            scale_minor,
            center,
            translation,
        }
    }
    pub fn reverse(mut self) -> Self {
        self.full_angle *= -1.0;
        self.major_axis = rotate(self.cos, self.sin, self.major_axis);
        self.sin *= -1.0;
        self.center = self.center + self.translation;
        self.translation *= -1.0;
        self.scale_major = 1.0 / self.scale_major;
        self.scale_minor = 1.0 / self.scale_minor;

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
        let minor_axis = self.major_axis.perp();
        let v = self.major_axis * self.major_axis.dot(rhs) * self.scale_major
            + minor_axis * minor_axis.dot(rhs) * self.scale_minor;
        self.center + self.translation + rotate(self.cos, self.sin, v)
    }
}

fn rotate(cos: f32, sin: f32, v: Vec2) -> Vec2 {
    Vec2::new(v.x * cos - v.y * sin, v.y * cos + v.x * sin)
}
