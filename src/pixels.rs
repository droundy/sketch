#[derive(Default, Clone, PartialEq, Eq)]
pub struct Pixels {
    ranges: Vec<Range>,
}

#[derive(Default, Clone, PartialEq, Eq)]
struct Range {
    start: u32,
    length: u32,
}

impl Pixels {
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.ranges
            .iter()
            .flat_map(|r| r.start as usize..(r.start + r.length) as usize)
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
                    ranges.push(Range {
                        start: start as u32,
                        length: (i - start) as u32,
                    });
                    started_at = None;
                }
            }
        }
        if let Some(start) = started_at {
            ranges.push(Range {
                start: start as u32,
                length: (bits.len() - start) as u32,
            })
        }
        Pixels { ranges }
    }
}
