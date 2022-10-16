use std::{ops::Mul, time::Instant};

use crate::pixeltree::Pixels;
use glam::Vec2;
use serde::{Deserialize, Serialize};

/// A connection between two keyframes.
#[derive(Clone, Serialize, Deserialize)]
pub struct Tween {
    chunks_alternating_fill: Vec<Vec<ChunkTween>>,
}

impl Tween {
    pub fn new(w: usize, before: Pixels, after: Pixels) -> Self {
        if before.is_empty() || after.is_empty() {
            return Tween {
                chunks_alternating_fill: Vec::new(),
            };
        }
        let before_copy = before.clone();
        let after_copy = after.clone();

        let mut before_chunks = Chunk::find(w, before);
        before_chunks.sort_by(|a, b| b.area.cmp(&a.area));
        let mut after_chunks = Chunk::find(w, after);
        after_chunks.sort_by(|a, b| b.area.cmp(&a.area));

        let mut before_matched = Vec::new();
        let mut after_matched = Vec::new();
        // First match up any chunks that exist unchanged before and after
        for b in (0..before_chunks.len()).rev() {
            for a in 0..after_chunks.len() {
                if before_chunks[b].points == after_chunks[a].points {
                    before_matched.push(before_chunks.remove(b));
                    after_matched.push(after_chunks.remove(a));
                    break;
                }
            }
        }
        // Then match up chunks that have only been moved but not changed.
        for b in (0..before_chunks.len()).rev() {
            if before_chunks[b + 1..]
                .iter()
                .filter(|c| c.shifted_eq(&before_chunks[b]))
                .next()
                .is_none()
                && after_chunks
                    .iter()
                    .filter(|c| c.shifted_eq(&before_chunks[b]))
                    .count()
                    == 1
            {
                for a in 0..after_chunks.len() {
                    if before_chunks[b].shifted_eq(&after_chunks[a]) {
                        before_matched.push(before_chunks.remove(b));
                        after_matched.push(after_chunks.remove(a));
                        break;
                    }
                }
            }
        }

        let nchunks = std::cmp::min(before_chunks.len(), after_chunks.len());
        while before_chunks.len() > nchunks {
            before_chunks.pop();
        }
        while after_chunks.len() > nchunks {
            after_chunks.pop();
        }
        if nchunks != 0 {
            // FIXME this is a hokwy way to pair up the chunks.
            for _ in 0..10000 {
                let i = rand::random::<usize>() % nchunks;
                let j = rand::random::<usize>() % nchunks;
                let v_before = before_chunks[i].center - before_chunks[j].center;
                let v_after = after_chunks[i].center - after_chunks[j].center;
                if v_after.dot(v_before) < 0.0 {
                    after_chunks.swap(i, j);
                }
            }
        }
        before_chunks.extend(before_matched);
        after_chunks.extend(after_matched);
        let mut chunks_alternating_fill = vec![Vec::new()];
        for (before, after) in before_chunks.into_iter().zip(after_chunks.into_iter()) {
            let mut before_fill: Pixels = before.points.compute_fill(w);
            before_fill.remove(&before_copy);
            let mut after_fill: Pixels = after.points.compute_fill(w);
            after_fill.remove(&after_copy);
            let main_tween = ChunkTween::new(w, before, after);
            chunks_alternating_fill[0].push(main_tween);
            if !before_fill.is_empty() && !after_fill.is_empty() {
                let tween_fill = Tween::new(w, before_fill, after_fill);
                // let tween_fill = Tween::with_transform(w, before_fill, after_fill);
                while chunks_alternating_fill.len() < tween_fill.chunks_alternating_fill.len() + 1 {
                    chunks_alternating_fill.push(Vec::new());
                }
                for (c, add) in chunks_alternating_fill
                    .iter_mut()
                    .skip(1)
                    .zip(tween_fill.chunks_alternating_fill.into_iter())
                {
                    c.extend(add);
                }
            }
        }
        Tween {
            chunks_alternating_fill,
        }
    }

    /// Draw this Tween to the pixel buffer
    pub fn draw(&self, fraction: f32, pixels: &mut [bool]) {
        for p in self.draw_points(fraction).iter() {
            if let Some(x) = pixels.get_mut(p) {
                *x = true;
            }
        }
    }

    /// Draw this Tween to the pixel buffer
    pub fn draw_points(&self, fraction: f32) -> Pixels {
        let mut final_pixels = Pixels::default();
        let w = if let Some(Some(c)) = self.chunks_alternating_fill.first().map(|x| x.first()) {
            c.w
        } else {
            return final_pixels;
        };
        for chunk_and_fill in self.chunks_alternating_fill.chunks(2) {
            let mut outline = Pixels::default();
            let mut fill = Pixels::default();
            for c in chunk_and_fill[0].iter() {
                outline.extend(&c.interpolate(fraction));
            }
            let mut points = outline.clone();
            points.extend(&outline.compute_fill(w));
            if let Some(filltweens) = chunk_and_fill.get(1) {
                // This means we are not the last odd chunk, and need to remove our fill.
                for c in filltweens.iter() {
                    fill.extend(&c.interpolate(fraction));
                }
                points.remove(&fill);
                points.remove(&fill.compute_fill(w));
                points.extend(&outline);
            }
            final_pixels.extend(&points)
        }
        final_pixels
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
        let mut after_pixels = after.points.clone();

        let before_outline = outline(w, &mut before_pixels).unwrap();
        let after_outline = outline(w, &mut after_pixels).unwrap();

        let before_positions = before_outline
            .iter()
            .map(|i| transform * Vec2::new((i % w) as f32, (i / w) as f32))
            .collect::<Vec<_>>();
        let after_positions = after_outline
            .iter()
            .map(|i| Vec2::new((i % w) as f32, (i / w) as f32))
            .collect::<Vec<_>>();

        let (smaller_positions, larger_positions) = if before_outline.len() < after_outline.len() {
            (before_positions, after_positions)
        } else {
            (after_positions, before_positions)
        };

        let mut outline_connections = Vec::with_capacity(smaller_positions.len());
        // We start by identifying the closest point on the outlines by distance.
        for (i, pos) in smaller_positions.iter().copied().enumerate() {
            let mut j = 0;
            let mut closest = larger_positions[j].distance_squared(pos);
            for (jj, posj) in larger_positions.iter().copied().enumerate() {
                if posj.distance_squared(pos) < closest {
                    j = jj;
                    closest = posj.distance_squared(pos);
                }
            }
            outline_connections.push((i, j));
        }
        // We then drop all the connections that are out of order.
        let outline_connections = make_monotonic(outline_connections);
        // Now reorder the indexes if needed
        let outline_connections = if before_outline.len() < after_outline.len() {
            outline_connections
        } else {
            outline_connections.iter().map(|&(i, j)| (j, i)).collect()
        };
        // Finally we interpolate the outline in between those closest connections, so the
        // whole outline should continuously deform with no breaks.
        let (mut i, mut j) = outline_connections.last().copied().unwrap();
        for (mut b, mut a) in outline_connections {
            if b < i {
                b += before_outline.len();
            }
            if a < j {
                a += after_outline.len();
            }
            let gap_b = b - i + 1;
            let gap_a = a - j + 1;
            let initial_i = i;
            let initial_j = j;
            connections.push((
                before_outline[i % before_outline.len()],
                after_outline[j % after_outline.len()],
            ));
            while i < b || j < a {
                // println!("     ({i}, {j}) {}", connections.len());
                connections.push((
                    before_outline[i % before_outline.len()],
                    after_outline[j % after_outline.len()],
                ));
                let ifrac = (i + 1 - initial_i) as f64 / gap_b as f64;
                let jfrac = (j + 1 - initial_j) as f64 / gap_a as f64;
                if ifrac < jfrac || j >= a {
                    i += 1;
                } else {
                    j += 1;
                }
            }
            i = i % before_outline.len();
            j = j % after_outline.len();
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
    fn interpolate(&self, fraction: f32) -> Pixels {
        assert!(fraction >= 0.0);
        assert!(fraction <= 1.0);
        let reverse_transform = (1.0 - fraction) * self.transform.reverse();
        let transform = fraction * self.transform;
        let start = Instant::now();
        let mut pixels = Pixels::default();
        for &(b, a) in self.connections.iter() {
            let b = transform * Vec2::new((b % self.w) as f32, (b / self.w) as f32);
            let a = reverse_transform * Vec2::new((a % self.w) as f32, (a / self.w) as f32);
            let p = fraction * a + (1.0 - fraction) * b;
            let w = self.w;
            let idx0 = p.x as usize + (p.y as usize) * w;
            pixels.insert(idx0);
            pixels.insert(idx0 + 1);
            pixels.insert(idx0 + w);
            pixels.insert(idx0 + w + 1);
        }
        if start.elapsed().as_secs_f64() > 6e-3 {
            println!(
                "Number connections is {}, took {:.2} ms",
                self.connections.len(),
                start.elapsed().as_secs_f64() * 1e3
            );
        }
        pixels
    }
}

#[derive(Clone)]
pub struct Chunk {
    points: Pixels,
    center: Vec2,
    extrema: [Vec2; 4],
    area: usize,
    // The major axis length
    major: f32,
    // The minor axis length
    minor: f32,
    // The direction of the major axis.
    axis: Vec2,
}

impl Chunk {
    pub fn find(w: usize, mut pixels: Pixels) -> Vec<Self> {
        let mut out = Vec::new();
        while let Some(points) = contiguous_pixels(w, &mut pixels) {
            out.push(Chunk::new(w, points));
        }
        out
    }
    pub fn new(w: usize, points: Pixels) -> Self {
        let mut center = Vec2::ZERO;
        let top = points.iter().next().unwrap();
        let mut top = Vec2::new((top % w) as f32, (top / w) as f32);
        let mut left = top;
        let mut bottom = top;
        let mut right = top;
        let mut area = 0;
        for p in points.iter() {
            area += 1;
            let v = Vec2::new((p % w) as f32, (p / w) as f32);
            center += v;
            // Combine x and y to get approximations of x and y that are likely
            // to be unique, so we will always get the same extrema, regardless
            // of the order of iteration.
            let skewy = |v: Vec2| v.y as f64 + 1e-7 * v.x as f64;
            let skewx = |v: Vec2| v.x as f64 + 1e-7 * v.y as f64;
            if skewy(v) < skewy(top) {
                top = v;
            }
            if skewy(v) > skewy(bottom) {
                bottom = v;
            }
            if skewx(v) < skewx(left) {
                left = v;
            }
            if skewx(v) > skewx(right) {
                right = v;
            }
        }
        let extrema = [top, right, bottom, left];
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
        x2 /= area as f32;
        y2 /= area as f32;
        xy /= area as f32;
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
            extrema,
            major,
            minor,
            axis,
        }
    }

    fn shifted_eq(&self, other: &Chunk) -> bool {
        self.points.shifted_eq(&other.points)
    }
}

fn contiguous_pixels(w: usize, pixels: &mut Pixels) -> Option<Pixels> {
    let mut out = Pixels::default();

    let p = pixels.iter().next()?;
    pixels.remove_pixel(p);
    let mut todo = vec![p];
    while let Some(p) = todo.pop() {
        out.insert(p);
        if p > 0 && pixels.contains(p - 1) {
            todo.push(p - 1);
            pixels.remove_pixel(p - 1);
        }
        if pixels.contains(p + 1) {
            todo.push(p + 1);
            pixels.remove_pixel(p + 1);
        }
        if p >= w && pixels.contains(p - w) {
            todo.push(p - w);
            pixels.remove_pixel(p - w);
        }
        if pixels.contains(p + w) {
            todo.push(p + w);
            pixels.remove_pixel(p + w);
        }
    }
    Some(out)
}

fn outline(w: usize, pixels: &mut Pixels) -> Option<Vec<usize>> {
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
                last
            }
        } else if next == last + w {
            if pixels.contains(next + 1) {
                next + 1
            } else if pixels.contains(next + w) {
                next + w
            } else if pixels.contains(next.wrapping_sub(1)) {
                next.wrapping_sub(1)
            } else {
                last
            }
        } else if next == last.wrapping_sub(1) {
            if pixels.contains(next + w) {
                next + w
            } else if pixels.contains(next.wrapping_sub(1)) {
                next.wrapping_sub(1)
            } else if pixels.contains(next.wrapping_sub(w)) {
                next.wrapping_sub(w)
            } else {
                last
            }
        } else if next == last.wrapping_sub(w) {
            if pixels.contains(next.wrapping_sub(1)) {
                next.wrapping_sub(1)
            } else if pixels.contains(next.wrapping_sub(w)) {
                next.wrapping_sub(w)
            } else if pixels.contains(next + 1) {
                next + 1
            } else {
                last
            }
        } else {
            unreachable!()
        };
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
    let pixels = Pixels::from_iter([2, 5, 8]);
    let w = 3;
    let chunks = Chunk::find(w, pixels);
    assert_eq!(1, chunks.len());
    assert_eq!(3, chunks[0].area);
    assert_eq!(Vec2::new(2.0, 1.0), chunks[0].center);
    assert_eq!(Vec2::new(0.0, 1.0), chunks[0].axis);
    assert_eq!(2.0_f32.sqrt(), chunks[0].major);
    assert_eq!(0.0, chunks[0].minor);

    let pixels = Pixels::from_iter([2, 3, 5, 6, 8]);
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
        let num_equal_extrema = o
            .extrema
            .iter()
            .zip(n.extrema.iter())
            .filter(|(a, b)| a == b)
            .count();

        let differences = |extrema: [Vec2; 4]| {
            [
                extrema[1] - extrema[0],
                extrema[2] - extrema[1],
                extrema[3] - extrema[2],
                extrema[0] - extrema[3],
            ]
        };
        let num_equal_differences = differences(o.extrema)
            .into_iter()
            .zip(differences(n.extrema).into_iter())
            .filter(|(a, b)| !a.abs_diff_eq(Vec2::ZERO, 0.1) && a.abs_diff_eq(*b, 0.5))
            .count();
        let mut full_angle = if num_equal_differences >= 1 {
            // If the vector difference between two extrema has not changed, then we can conclude that
            // the chunk as a whole has not rotated.
            0.0
        } else if o.axis.dot(n.axis) > 0.0 {
            -n.axis.angle_between(o.axis)
        } else {
            -n.axis.angle_between(-o.axis)
        };
        // println!("full_angle is {full_angle}");
        let (sin, cos) = full_angle.sin_cos();
        let major_axis = o.axis;
        let mut scale_major = if full_angle == 0.0 && o.axis.angle_between(n.axis).abs() > 0.05 {
            // Do not scale the major or minor axis, if we have eliminated rotation, and yet
            // the major/minor axes have shifted.
            1.0
        } else if o.major > 0.0 && n.major > 0.0 {
            n.major / o.major
        } else {
            1.0
        };
        let mut scale_minor = if full_angle == 0.0 && o.axis.angle_between(n.axis).abs() > 0.05 {
            // Do not scale the major or minor axis, if we have eliminated rotation, and yet
            // the major/minor axes have shifted.
            1.0
        } else if o.minor > 0.0 && n.minor > 0.0 {
            n.minor / o.minor
        } else {
            1.0
        };
        if o.major / o.minor < 1.1 || n.major / n.minor < 1.1 {
            if full_angle.abs() > 0.1 {
                let scale = (scale_major * scale_minor).sqrt();
                scale_major = scale;
                scale_minor = scale;
            }
            full_angle = 0.0;
        }
        // println!("scales are: {scale_major} and {scale_minor}");
        // println!("old major {} and minor {}", o.major, o.minor);
        // println!("new major {} and minor {}", n.major, n.minor);
        // println!("the major axis is {} {}", major_axis.x, major_axis.y);
        let mut center = o.center;
        let mut translation = n.center - o.center;
        if num_equal_extrema > 0 {
            translation = Vec2::ZERO;
        } else if full_angle == 0.0 && num_equal_differences > 0 {
            // We want to set the center to be the center of the rigid part
            // of the image.

            let which_pair_equal = differences(o.extrema)
                .into_iter()
                .zip(differences(n.extrema).into_iter())
                .enumerate()
                .filter(|(_, (a, b))| !a.abs_diff_eq(Vec2::ZERO, 0.1) && a.abs_diff_eq(*b, 0.5))
                .map(|(i, _)| i)
                .next()
                .unwrap();
            translation = n.extrema[which_pair_equal] - o.extrema[which_pair_equal];
            center = o.extrema[which_pair_equal];
        }
        // println!("full_angle {full_angle}");
        // println!("num_equal_extrema {num_equal_extrema}");
        // println!("num_equal_differences {num_equal_differences}");
        // for i in 0..4 {
        //     println!("    {:?} {:?}", o.extrema[i], n.extrema[i]);
        // }
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

fn is_monotonic(connections: &[(usize, usize)]) -> bool {
    count_monotonic(
        connections
            .iter()
            .copied()
            .chain(connections[0..1].iter().copied()),
        connections[0].1,
    ) == connections.len() + 1
}

#[test]
fn test_is_monotonic() {
    assert!(is_monotonic(&[(0, 1), (0, 2), (0, 4)]));
    assert!(is_monotonic(&[(0, 1), (0, 2), (0, 4), (0, 0)]));
    assert!(is_monotonic(&[(0, 1), (0, 2), (0, 3), (0, 0)]));
    assert!(is_monotonic(&[(0, 3), (0, 4), (0, 1), (0, 2)]));
    assert!(!is_monotonic(&[(0, 1), (0, 0), (0, 4)]));
    assert!(!is_monotonic(&[(0, 1), (0, 2), (0, 3), (0, 1), (0, 0)]));
}

fn count_monotonic(mut connections: impl Iterator<Item = (usize, usize)>, first: usize) -> usize {
    let mut prev = if let Some(prev) = connections.next() {
        prev.1
    } else {
        return 0;
    };
    let mut have_wrapped = prev < first;
    let mut count = 1;
    for (_, next) in connections {
        if have_wrapped {
            if next < prev || next > first {
                return count;
            }
        } else {
            if next < prev {
                if next > first {
                    return count;
                }
                have_wrapped = true;
            }
        }
        prev = next;
        count += 1;
    }
    count
}

#[test]
fn test_count_monotonic() {
    assert_eq!(
        3,
        count_monotonic([(0, 1), (0, 0), (0, 1), (0, 2), (0, 3)].into_iter(), 1)
    );
    assert_eq!(
        5,
        count_monotonic([(0, 1), (0, 2), (0, 2), (0, 2), (0, 3)].into_iter(), 1)
    );
}

fn make_monotonic(mut connections: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    let mut is_keeper = vec![false; connections.len()];
    while !is_monotonic(&connections) {
        let mut longest_run = 0..0;
        for i in 0..connections.len() {
            if !is_keeper[i] {
                let run = count_monotonic(
                    connections[i..]
                        .iter()
                        .copied()
                        .chain(connections[0..i + 1].iter().copied()),
                    connections[i].1,
                );
                if run > longest_run.len() {
                    longest_run = i..i + run;
                }
            }
        }
        if longest_run.len() < 2 {
            // FIXME
            break;
        }
        for i in longest_run.clone() {
            is_keeper[i % connections.len()] = true;
        }
        // Now we need to eliminate any connections that have already been ruled out.
        // how we do this will depend on whether or not we have wrapped around the
        // end in our monotonic run.
        let run_end_value = connections[longest_run.end % connections.len()].1;
        let run_start_value = connections[longest_run.start].1;
        let max_to_keep = std::cmp::min(run_end_value, run_start_value);
        let min_to_keep = std::cmp::max(run_end_value, run_start_value);

        let mut to_remove = vec![false; connections.len()];
        for j in longest_run.end..longest_run.end + connections.len() {
            let j = j % connections.len();
            if is_keeper[j] {
                break;
            }
            let v = connections[j].1;
            if v > max_to_keep || v < min_to_keep {
                to_remove[j] = true;
            }
        }
        for (i, _) in to_remove
            .iter()
            .copied()
            .enumerate()
            .rev()
            .filter(|(_, b)| *b)
        {
            is_keeper.remove(i);
            connections.remove(i);
        }
    }
    assert!(is_monotonic(&connections));
    connections
}

#[test]
fn test_make_monotonic() {
    assert_eq!(vec![(0, 1), (0, 2)], make_monotonic(vec![(0, 1), (0, 2)]));
    assert_eq!(
        vec![(0, 1), (0, 2), (0, 3), (0, 0)],
        make_monotonic(vec![(0, 1), (0, 2), (0, 3), (0, 0)])
    );
    assert_eq!(
        vec![(0, 1), (0, 2), (0, 3), (0, 1)],
        make_monotonic(vec![(0, 1), (0, 2), (0, 3), (0, 1), (0, 0)])
    );
}
