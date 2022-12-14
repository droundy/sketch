use std::{borrow::Borrow, collections::BTreeSet, iter::FromIterator};

use serde::{Deserialize, Serialize};

#[derive(Default, Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Pixels {
    ranges: BTreeSet<Borders>,
}

/// These ranges are equal if they overlap or border one another.
#[derive(Default, Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
struct Borders(Overlaps);

impl Borrow<Overlaps> for Borders {
    fn borrow(&self) -> &Overlaps {
        &self.0
    }
}

impl PartialOrd for Borders {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.0.start + self.0.length + 1 < other.0.start {
            Some(std::cmp::Ordering::Less)
        } else if other.0.start + other.0.length + 1 < self.0.start {
            Some(std::cmp::Ordering::Greater)
        } else {
            Some(std::cmp::Ordering::Equal)
        }
    }
}
impl Ord for Borders {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.0.start + self.0.length < other.0.start {
            std::cmp::Ordering::Less
        } else if other.0.start + other.0.length < self.0.start {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

/// These ranges are equal if one overlaps the other.
#[derive(Default, Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
struct Overlaps {
    start: i32,
    length: i32,
}

impl PartialOrd for Overlaps {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.start + self.length - 1 < other.start {
            Some(std::cmp::Ordering::Less)
        } else if other.start + other.length - 1 < self.start {
            Some(std::cmp::Ordering::Greater)
        } else {
            Some(std::cmp::Ordering::Equal)
        }
    }
}
impl Ord for Overlaps {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.start + self.length - 1 < other.start {
            std::cmp::Ordering::Less
        } else if other.start + other.length - 1 < self.start {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

impl Pixels {
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.ranges
            .iter()
            .filter(|r| r.0.start + r.0.length >= 0)
            .flat_map(|r| r.0.start..(r.0.start + r.0.length))
            .flat_map(|p| p.try_into())
    }
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }
    pub fn contains(&self, elem: usize) -> bool {
        self.ranges.contains(&Overlaps {
            start: elem as i32,
            length: 1,
        })
    }
    pub fn contiguous_with(&self, w: usize, other: &Pixels) -> Pixels {
        let w = w as i32;
        let mut temp = self.clone();
        let mut contiguous = Pixels::default();
        let mut todo = Vec::new();
        for r in other.ranges.iter().map(|b| b.0) {
            while let Some(b) = temp.ranges.take(&r) {
                contiguous.ranges.insert(b);
                todo.push(b);
            }
        }
        while let Some(r) = todo.pop() {
            while let Some(b) = temp.ranges.take(&r) {
                contiguous.ranges.insert(b);
                todo.push(b);
            }
            if r.0.start > w {
                let up = Overlaps {
                    start: r.0.start - w,
                    length: r.0.length,
                };
                while let Some(b) = temp.ranges.take(&up) {
                    contiguous.ranges.insert(b);
                    todo.push(b);
                }
            }
            let down = Overlaps {
                start: r.0.start + w,
                length: r.0.length,
            };
            while let Some(b) = temp.ranges.take(&down) {
                contiguous.ranges.insert(b);
                todo.push(b);
            }
        }
        contiguous
    }
    pub fn overlaps(&self, other: &Pixels) -> bool {
        if other.ranges.len() < self.ranges.len() {
            other.overlaps(self)
        } else {
            for r in self.ranges.iter().map(|x| &x.0) {
                if other.ranges.contains(r) {
                    return true;
                }
            }
            false
        }
    }
    pub fn borders(&self, other: &Pixels) -> bool {
        if other.ranges.len() < self.ranges.len() {
            other.borders(self)
        } else {
            for r in self.ranges.iter() {
                if other.ranges.contains(r) {
                    return true;
                }
            }
            false
        }
    }
    pub fn shift_by(&mut self, offset: isize) {
        let offset = offset as i32;
        let old_ranges = std::mem::take(&mut self.ranges);
        self.ranges = old_ranges
            .into_iter()
            .map(|mut r| {
                r.0.start += offset;
                r
            })
            .collect();
    }
    pub fn shifted_eq(&self, other: &Pixels) -> bool {
        if self.ranges.len() != other.ranges.len() {
            return false;
        }
        if self.ranges.is_empty() {
            return true;
        }
        let offset = self.ranges.iter().next().unwrap().0.start
            - other.ranges.iter().next().unwrap().0.start;
        for (rself, rother) in self.ranges.iter().zip(other.ranges.iter()) {
            if rself.0.length != rother.0.length || rself.0.start != rother.0.start + offset {
                return false;
            }
        }
        true
    }
    pub fn insert(&mut self, elem: usize) {
        self.insert_borders(Borders(Overlaps {
            start: elem as i32,
            length: 1,
        }))
    }
    fn insert_borders(&mut self, mut b: Borders) {
        while !self.ranges.insert(b) {
            let o = self.ranges.take(&b).unwrap();
            let start = std::cmp::min(b.0.start, o.0.start);
            let stop = std::cmp::max(o.0.start + o.0.length, b.0.start + b.0.length);
            let length = stop - start;
            b = Borders(Overlaps { start, length });
        }
    }
    pub fn extend(&mut self, pixels: &Pixels) {
        for b in pixels.ranges.iter().copied() {
            self.insert_borders(b);
        }
    }
    pub fn expand(&self, w: usize) -> Pixels {
        let mut out = Pixels::default();
        let w = w as i32;
        for r in self.ranges.iter().copied() {
            out.insert_borders(Borders(Overlaps {
                start: r.0.start - w,
                length: r.0.length,
            }));
            out.insert_borders(Borders(Overlaps {
                start: r.0.start - 1,
                length: r.0.length + 2,
            }));
            out.insert_borders(Borders(Overlaps {
                start: r.0.start + w,
                length: r.0.length,
            }));
        }
        out
    }
    pub fn expansion_region(&self, w: usize) -> Pixels {
        let mut out = Pixels::default();
        let w = w as i32;
        for r in self.ranges.iter().copied() {
            out.insert_borders(Borders(Overlaps {
                start: r.0.start - w,
                length: r.0.length,
            }));
            out.insert_borders(Borders(Overlaps {
                start: r.0.start - 1,
                length: 1,
            }));
            out.insert_borders(Borders(Overlaps {
                start: r.0.start + r.0.length,
                length: 1,
            }));
            out.insert_borders(Borders(Overlaps {
                start: r.0.start + w,
                length: r.0.length,
            }));
        }
        out.remove(self);
        out
    }
    fn remove_borders(&mut self, b: Borders) {
        while let Some(o) = self.ranges.take(&b.0) {
            // No need for the more expensive insert_borders, since we know there
            // will not be an overlap.
            if o.0.start < b.0.start {
                self.ranges.insert(Borders(Overlaps {
                    length: b.0.start - o.0.start,
                    ..o.0
                }));
            }
            if o.0.start + o.0.length > b.0.start + b.0.length {
                self.ranges.insert(Borders(Overlaps {
                    length: o.0.start + o.0.length - b.0.start - b.0.length,
                    start: b.0.start + b.0.length,
                }));
            }
        }
    }

    pub fn remove(&mut self, pix: &Pixels) {
        if pix.ranges.len() >= self.ranges.len() {
            for b in pix.ranges.iter().copied() {
                self.remove_borders(b);
            }
        } else {
            let mut to_remove = Vec::new();
            let mut found_overlap = true;
            while found_overlap {
                for b in self.ranges.iter() {
                    found_overlap = false;
                    if let Some(bb) = pix.ranges.get(&b.0).copied() {
                        to_remove.push(bb);
                        found_overlap = true;
                    }
                }
                for b in to_remove.drain(..) {
                    self.remove_borders(b);
                }
            }
        }
    }
    pub fn remove_pixel(&mut self, p: usize) {
        self.remove_borders(Borders(Overlaps {
            start: p as i32,
            length: 1,
        }));
    }

    pub fn keep_intersection(&mut self, pix: &Pixels) {
        let mut mine = std::mem::take(&mut self.ranges);
        for r in pix.ranges.iter().map(|b| b.0) {
            while let Some(rr) = mine.take(&r) {
                let start = std::cmp::max(rr.0.start, r.start);
                let end = std::cmp::min(rr.0.start + rr.0.length, r.start + r.length);
                let length = end - start;
                self.ranges.insert(Borders(Overlaps { start, length }));
            }
        }
    }

    pub fn compute_fill(&self, w: usize) -> Pixels {
        let w = w as i32;
        let min = if let Some(min) = self.ranges.iter().next().copied() {
            min.0.start
        } else {
            return Pixels::default();
        };
        let max_plus_1 = if let Some(max) = self.ranges.iter().rev().next().copied() {
            max.0.start + max.0.length
        } else {
            return Pixels::default();
        };
        let mut inside = Pixels {
            ranges: [Borders(Overlaps {
                start: min - w,
                length: max_plus_1 - min + 2 * w,
            })]
            .into_iter()
            .collect(),
        };
        inside.remove(self);
        let mut todo = vec![
            Overlaps {
                start: min - w - 1,
                length: w,
            },
            Overlaps {
                start: max_plus_1 + w,
                length: w,
            },
        ];
        while let Some(outside) = todo.pop() {
            while let Some(out) = inside.ranges.take(&outside) {
                todo.push(Overlaps {
                    start: out.0.start + w,
                    length: out.0.length,
                });
                todo.push(Overlaps {
                    start: out.0.start - w,
                    length: out.0.length,
                });
            }
        }
        inside
    }

    pub fn count(&self) -> usize {
        self.ranges.iter().map(|r| r.0.length as usize).sum()
    }
}

impl FromIterator<usize> for Pixels {
    fn from_iter<T: IntoIterator<Item = usize>>(iter: T) -> Self {
        let mut iter = iter.into_iter();
        let mut out = Pixels {
            ranges: std::collections::BTreeSet::new(),
        };
        if let Some(mut prev) = iter.next() {
            let mut r = Borders(Overlaps {
                start: prev as i32,
                length: 1,
            });
            for p in iter {
                if p == prev + 1 {
                    r.0.length += 1;
                } else {
                    out.insert_borders(r);
                    r = Borders(Overlaps {
                        start: p as i32,
                        length: 1,
                    });
                }
                prev = p;
            }
            out.insert_borders(r);
        }
        out
    }
}

impl From<Vec<bool>> for Pixels {
    fn from(bits: Vec<bool>) -> Self {
        assert!(bits.len() < i32::MAX as usize);
        let mut ranges = Vec::new();
        let mut started_at = None;
        for (i, b) in bits.iter().copied().enumerate() {
            if b && started_at.is_none() {
                started_at = Some(i);
            }
            if !b {
                if let Some(start) = started_at {
                    ranges.push(Borders(Overlaps {
                        start: start as i32,
                        length: (i - start) as i32,
                    }));
                    started_at = None;
                }
            }
        }
        if let Some(start) = started_at {
            ranges.push(Borders(Overlaps {
                start: start as i32,
                length: (bits.len() - start) as i32,
            }));
        }
        let ranges = ranges.into_iter().collect();
        Pixels { ranges }
    }
}

#[test]
fn testme() {
    let mut v1 = vec![false, false, true, true, false, true];
    let mut set1 = Pixels::from(v1.clone());
    println!("set1 is {set1:?}");
    for i in 0..v1.len() {
        println!("Element {i} should be {}", v1[i]);
        assert_eq!(v1[i], set1.contains(i));
    }

    let set2 = Pixels::from(vec![
        false, false, false, false, false, false, true, true, false, true,
    ]);
    assert!(!set2.overlaps(&set1));
    assert!(!set1.overlaps(&set2));

    let set2 = Pixels::from(vec![false, false, false, true, true, false]);
    assert!(set2.overlaps(&set1));
    assert!(set1.overlaps(&set2));

    set1.insert(2);
    println!("set1 after inserting 2 is {set1:?}");
    for i in 0..v1.len() {
        println!("Element {i} should be {}", v1[i]);
        assert_eq!(v1[i], set1.contains(i));
    }

    set1.insert(0);
    set1.insert(1);
    v1[0] = true;
    v1[1] = true;
    println!("set1 after inserting 0 is {set1:?}");
    for i in 0..v1.len() {
        println!("Element {i} should be {}", v1[i]);
        assert_eq!(v1[i], set1.contains(i));
    }

    let s1: Pixels = [1, 2, 5, 7].into_iter().collect();
    let mut s2: Pixels = (1..8).collect();
    println!("\nFresh start! with {s1:?} and {s2:?}\n");
    assert!(s1.overlaps(&s2));
    let s2_minus_s1: Pixels = [3, 4, 6].into_iter().collect();

    s2.remove(&s1);
    println!("s2 is now {s2:?}");
    println!("s2 should be {s2_minus_s1:?}");
    for i in 0..16 {
        assert_eq!(s2.contains(i), s2_minus_s1.contains(i));
    }

    for v0 in [false, true] {
        for v1 in [false, true] {
            for v2 in [false, true] {
                for v3 in [false, true] {
                    let v = vec![v0, v1, v2, v3];
                    let s1 = Pixels::from(v.clone());
                    for i in 0..4 {
                        assert_eq!(v[i], s1.contains(i));
                    }
                    for i in 0..4 {
                        let mut v = v.clone();
                        let mut s1 = s1.clone();
                        v[i] = false;
                        s1.remove_pixel(i);
                        for j in 0..4 {
                            assert_eq!(v[j], s1.contains(j));
                        }
                    }
                }
            }
        }
    }
}
