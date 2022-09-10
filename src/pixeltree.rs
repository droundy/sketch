use std::{borrow::Borrow, collections::BTreeSet, iter::FromIterator};

#[derive(Default, Clone, PartialEq, Eq, Debug)]
pub struct Pixels {
    ranges: BTreeSet<Borders>,
}

/// These ranges are equal if they overlap or border one another.
#[derive(Default, Clone, Copy, PartialEq, Eq, Debug)]
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
#[derive(Default, Clone, Copy, PartialEq, Eq, Debug)]
struct Overlaps {
    start: u32,
    length: u32,
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
            .flat_map(|r| r.0.start as usize..(r.0.start + r.0.length) as usize)
    }
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }
    pub fn contains(&self, elem: usize) -> bool {
        self.ranges.contains(&Overlaps {
            start: elem as u32,
            length: 1,
        })
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
    pub fn insert(&mut self, elem: usize) {
        self.insert_borders(Borders(Overlaps {
            start: elem as u32,
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
    pub fn expand(&self, w: usize) -> Pixels {
        let mut out = Pixels::default();
        let w = w as u32;
        for r in self.ranges.iter().copied() {
            if r.0.start > w {
                out.insert_borders(Borders(Overlaps {
                    start: r.0.start - w,
                    length: r.0.length,
                }));
            }
            if r.0.start > 0 {
                out.insert_borders(Borders(Overlaps {
                    start: r.0.start - 1,
                    length: r.0.length + 2,
                }));
            }
            out.insert_borders(Borders(Overlaps {
                start: r.0.start + w,
                length: r.0.length,
            }));
        }
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
        for b in pix.ranges.iter().copied() {
            self.remove_borders(b);
        }
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
                start: prev as u32,
                length: 1,
            });
            for p in iter {
                if p == prev + 1 {
                    r.0.length += 1;
                } else {
                    out.insert_borders(r);
                    r = Borders(Overlaps {
                        start: p as u32,
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
        assert!(bits.len() < u32::MAX as usize);
        let mut ranges = Vec::new();
        let mut started_at = None;
        for (i, b) in bits.iter().copied().enumerate() {
            if b && started_at.is_none() {
                started_at = Some(i);
            }
            if !b {
                if let Some(start) = started_at {
                    ranges.push(Borders(Overlaps {
                        start: start as u32,
                        length: (i - start) as u32,
                    }));
                    started_at = None;
                }
            }
        }
        if let Some(start) = started_at {
            ranges.push(Borders(Overlaps {
                start: start as u32,
                length: (bits.len() - start) as u32,
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
}
