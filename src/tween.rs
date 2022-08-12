use std::{
    collections::{HashMap, HashSet, VecDeque},
    ops::Mul,
};

use macroquad::prelude::Vec2;
use ordered_float::OrderedFloat;

/// A connection between two keyframes.
pub struct Tween {
    chunks: Vec<ChunkTween>,
}

impl Tween {
    pub fn new(w: usize, before: Vec<bool>, after: Vec<bool>) -> Self {
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
    ///
    /// FIXME: Need to respect the transform!
    pub fn draw(&self, fraction: f32, color: [u8; 4], pixels: &mut [[u8; 4]]) {
        for c in self.chunks.iter() {
            c.draw(fraction, color, pixels);
        }
    }
}

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

        let before_positions = before
            .points
            .iter()
            .copied()
            .map(|i| (i, Vec2::new((i % w) as f32, (i / w) as f32)))
            .map(|(i, v)| (i, transform * v))
            .collect::<Vec<_>>();
        let after_positions = after
            .points
            .iter()
            .copied()
            .map(|i| (i, Vec2::new((i % w) as f32, (i / w) as f32)))
            .collect::<Vec<_>>();

        let top_right_before = before_positions
            .iter()
            .max_by(|a, b| (a.1.x - a.1.y).partial_cmp(&(b.1.x - b.1.y)).unwrap())
            .unwrap()
            .0;
        let top_left_before = before_positions
            .iter()
            .max_by(|a, b| (-a.1.x - a.1.y).partial_cmp(&(-b.1.x - b.1.y)).unwrap())
            .unwrap()
            .0;

        let top_right_after = after_positions
            .iter()
            .max_by(|a, b| (a.1.x - a.1.y).partial_cmp(&(b.1.x - b.1.y)).unwrap())
            .unwrap()
            .0;
        let top_left_after = after_positions
            .iter()
            .max_by(|a, b| (-a.1.x - a.1.y).partial_cmp(&(-b.1.x - b.1.y)).unwrap())
            .unwrap()
            .0;

        let mut connections = Vec::new();

        let before_top_rankings = rank_pixels(w, before.to_pixels(), top_right_before);
        let before_left_rankings = rank_pixels(w, before.to_pixels(), top_left_before);

        assert_eq!(before_top_rankings.len(), before_left_rankings.len());

        let after_top_rankings = rank_pixels(w, after.to_pixels(), top_right_after);
        let after_left_rankings = rank_pixels(w, after.to_pixels(), top_left_after);

        assert_eq!(
            after_top_rankings.keys().copied().collect::<HashSet<_>>(),
            after_left_rankings.keys().copied().collect::<HashSet<_>>()
        );
        assert_eq!(
            before_top_rankings.keys().copied().collect::<HashSet<_>>(),
            before_left_rankings.keys().copied().collect::<HashSet<_>>()
        );

        let mut before_rankings = Vec::new();
        for (&k, &top) in before_top_rankings.iter() {
            let left = before_left_rankings[&k];
            before_rankings.push((top, left, k));
        }
        before_rankings.sort_unstable();

        let mut after_rankings = Vec::new();
        for (&k, &top) in after_top_rankings.iter() {
            let left = after_left_rankings[&k];
            after_rankings.push((top, left, k));
        }
        after_rankings.sort_unstable();

        for &(t, l, b) in before_rankings.iter() {
            let i = after_rankings.binary_search(&(t, l, b));
            let i = match i {
                Ok(i) => i,
                Err(i) => {
                    if i < after_rankings.len() && i > 0 {
                        let d_i =
                            (t - after_rankings[i].0).powi(2) + (l - after_rankings[i].1).powi(2);
                        let d_i_1 = (t - after_rankings[i - 1].0).powi(2)
                            + (l - after_rankings[i - 1].1).powi(2);
                        if d_i < d_i_1 {
                            i
                        } else {
                            i - 1
                        }
                    } else {
                        std::cmp::min(i, after_rankings.len() - 1)
                    }
                }
            };
            let (_, _, a) = after_rankings[i];
            connections.push((b, a))
        }
        for &(t, l, a) in after_rankings.iter() {
            let i = before_rankings.binary_search(&(t, l, a));
            let i = match i {
                Ok(i) => i,
                Err(i) => {
                    if i < before_rankings.len() && i > 0 {
                        let d_i =
                            (t - before_rankings[i].0).powi(2) + (l - before_rankings[i].1).powi(2);
                        let d_i_1 = (t - before_rankings[i - 1].0).powi(2)
                            + (l - before_rankings[i - 1].1).powi(2);
                        if d_i < d_i_1 {
                            i
                        } else {
                            i - 1
                        }
                    } else {
                        std::cmp::min(i, before_rankings.len() - 1)
                    }
                }
            };
            let (_, _, b) = before_rankings[i];
            connections.push((b, a))
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
    ///
    /// FIXME: Need to respect the transform!
    fn draw(&self, fraction: f32, color: [u8; 4], pixels: &mut [[u8; 4]]) {
        assert!(fraction >= 0.0);
        assert!(fraction <= 1.0);
        let reverse_transform = (1.0 - fraction) * self.transform.reverse();
        let transform = fraction * self.transform;
        for &(b, a) in self.connections.iter() {
            let b = transform * Vec2::new((b % self.w) as f32, (b / self.w) as f32);
            let a = reverse_transform * Vec2::new((a % self.w) as f32, (a / self.w) as f32);
            let p = fraction * a + (1.0 - fraction) * b;
            let idx = a.x.round() as usize + (a.y.round() as usize) * self.w;
            pixels[idx] = color;
            let idx = b.x.round() as usize - 100 + (b.y.round() as usize) * self.w;
            if idx < pixels.len() {
                pixels[idx] = [255, 0, 0, 255];
            }
        }
    }
}

pub struct Chunk {
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
            xy += -dx * dy;
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
            (-(y2 - x2) - ((x2 - y2).powi(2) + 4.0 * xy.powi(2)).sqrt()) / (2.0 * xy),
            1.0,
        )
        .normalize();
        let ax2 = Vec2::new(
            -((y2 - x2) - ((x2 - y2).powi(2) + 4.0 * xy.powi(2)).sqrt()) / (2.0 * xy),
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
        println!("[{x2:8.0} {xy:8.0}]   {:8.4e}   {:8.4e}", ax1.x, ax2.x);
        println!("[{xy:8.0} {y2:8.0}]   {:8.4e}   {:8.4e}", ax1.y, ax2.y);
        println!(
            "a1 mul moment is {} xhat + {} yhat",
            (x2 * ax1.x + xy * ax1.y) / ax1.x,
            (xy * ax1.x + y2 * ax1.y) / ax1.y
        );
        println!(
            "a2 mul moment is {} xhat + {} yhat",
            (x2 * ax2.x + xy * ax2.y) / ax2.x,
            (xy * ax2.x + y2 * ax2.y) / ax2.y
        );
        println!("{e1} {e2}");
        let (major, minor, axis) = if xy == 0.0 {
            if x2 > y2 {
                (x2.sqrt(), y2.sqrt(), Vec2::new(1.0, 0.0))
            } else {
                (y2.sqrt(), x2.sqrt(), Vec2::new(0.0, 1.0))
            }
        } else if e1 > e2 {
            (e1.sqrt(), e2.abs().sqrt(), ax1)
        } else {
            (e2.sqrt(), e1.abs().sqrt(), ax2)
        };
        Chunk {
            area,
            points,
            center,
            major,
            minor,
            axis,
        }
    }

    fn to_pixels(&self) -> Vec<bool> {
        let mut pixels = Vec::new();
        for p in self.points.iter().copied() {
            if p >= pixels.len() {
                pixels.extend(vec![false; p + 1]);
            }
            pixels[p] = true;
        }
        pixels
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
        println!("full_angle is {full_angle}");
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
        println!("scales are: {scale_major} and {scale_minor}");
        println!("old major {} and minor {}", o.major, o.minor);
        println!("new major {} and minor {}", n.major, n.minor);
        let center = o.center;
        let translation = n.center - o.center;
        Transform {
            full_angle,
            sin,
            cos,
            scale_major,
            major_axis,
            scale_minor,
            minor_axis,
            center,
            translation,
        }
    }
    pub fn reverse(mut self) -> Self {
        self.full_angle *= -1.0;
        self.major_axis = rotate(self.cos, self.sin, self.major_axis);
        self.minor_axis = rotate(self.cos, self.sin, self.minor_axis);
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
        let v = self.major_axis * self.major_axis.dot(rhs) * self.scale_major
            + self.minor_axis * self.minor_axis.dot(rhs) * self.scale_minor;
        self.center + self.translation + rotate(self.cos, self.sin, v)
    }
}

fn rotate(cos: f32, sin: f32, v: Vec2) -> Vec2 {
    Vec2::new(v.x * cos - v.y * sin, v.y * cos + v.x * sin)
}

fn rank_pixels(w: usize, mut pixels: Vec<bool>, start: usize) -> HashMap<usize, OrderedFloat<f64>> {
    assert!(pixels[start]);
    pixels[start] = false;
    let mut todo: VecDeque<(usize, i64)> = [(start, 0)].into_iter().collect();
    let mut out = HashMap::new();
    let mut max_rank = 0;
    while let Some((p, rank)) = todo.pop_front() {
        out.insert(p, rank);
        max_rank = rank;
        let x = p % w;
        let y = p / w;
        if x > 0 && pixels[p - 1] {
            todo.push_back((p - 1, rank + 1));
            pixels[p - 1] = false;
        }
        if x < w - 1 && pixels[p + 1] {
            todo.push_back((p + 1, rank + 1));
            pixels[p + 1] = false;
        }
        if y > 0 && pixels[p - w] {
            todo.push_back((p - w, rank + 1));
            pixels[p - w] = false;
        }
        if y < w - 1 && pixels[p + w] {
            todo.push_back((p + w, rank + 1));
            pixels[p + w] = false;
        }
    }
    let max_rank = max_rank as f64;
    out.into_iter()
        .map(|(k, v)| (k, (v as f64 / max_rank).into()))
        .collect()
}
