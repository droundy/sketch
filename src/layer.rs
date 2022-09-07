use std::collections::HashMap;

use macroquad::prelude::{screen_height, screen_width, Color, Image, Vec2};
use serde::{Deserialize, Serialize};
use tinyset::SetUsize;

use crate::{pixels::Pixels, tween::Tween};

#[derive(Clone, Serialize, Deserialize)]
struct Bitmap {
    time: f32,
    pixels: SetUsize,
    #[serde(skip)]
    fill_pixels: Pixels,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Layer {
    pub color: [u8; 4],
    pub fill_color: [u8; 4],
    width: usize,
    height: usize,
    // image: Image,
    // texture: Texture2D,
    keyframes: Vec<Bitmap>,
    #[serde(skip)]
    tweens: HashMap<(usize, usize), Tween>,
}

impl Layer {
    pub fn new(time: f32) -> Self {
        let bitmap = Image::gen_image_color(
            screen_width() as u16,
            screen_height() as u16,
            Color {
                r: 0.0,
                g: 0.0,
                b: 0.0,
                a: 0.0,
            },
        );
        Layer {
            color: random_color(),
            fill_color: [0; 4],
            width: bitmap.width(),
            height: bitmap.height(),
            // image: bitmap.clone(),
            // texture: Texture2D::from_image(&bitmap),
            keyframes: vec![Bitmap {
                time,
                // texture: Texture2D::from_image(&bitmap),
                pixels: SetUsize::new(),
                fill_pixels: Pixels::default(),
            }],
            tweens: HashMap::new(),
        }
    }
    pub fn get_color(&self) -> Color {
        Color::from_rgba(self.color[0], self.color[1], self.color[2], self.color[3])
    }
    pub fn get_fill_color(&self) -> Color {
        Color::from_rgba(
            self.fill_color[0],
            self.fill_color[1],
            self.fill_color[2],
            self.fill_color[3],
        )
    }
    fn compute_fill(&mut self, which: usize) {
        self.tweens.retain(|k, _| k.0 != which && k.1 != which);
        let k = &mut self.keyframes[which];
        {
            let mut inside = vec![true; self.width * self.height];
            let image_len = inside.len();
            for i in k.pixels.iter().filter(|&i| i < image_len) {
                inside[i] = false;
            }

            let mut todo = vec![0];
            let w = self.width;
            inside[0] = false;
            while let Some(i) = todo.pop() {
                if i > 0 && inside[i - 1] {
                    inside[i - 1] = false;
                    todo.push(i - 1);
                }
                if i + 1 < image_len && inside[i + 1] {
                    inside[i + 1] = false;
                    todo.push(i + 1);
                }
                if i >= w && inside[i - w] {
                    inside[i - w] = false;
                    todo.push(i - w);
                }
                if i + w < image_len && inside[i + w] {
                    inside[i + w] = false;
                    todo.push(i + w);
                }
            }
            k.fill_pixels = Pixels::from(inside);
        }
    }
    pub fn handle_modified_bitmap(&mut self, time: f32) {
        let which = self.closest_frame(time);
        self.compute_fill(which);
        let k = &mut self.keyframes[which];
        let mut img = Image::gen_image_color(
            self.width as u16,
            self.height as u16,
            Color::from_rgba(0, 0, 0, 0),
        );
        {
            let img = img.get_image_data_mut();
            let image_len = img.len();
            for i in k.pixels.iter().filter(|&i| i < image_len) {
                img[i] = self.color;
            }
            if self.fill_color[3] > 0 {
                for p in k.fill_pixels.iter() {
                    img[p] = self.fill_color;
                }
            }
        }
    }
    pub fn draw(&mut self, time: f32, pixels: &mut [[u8; 4]]) {
        let (before, after) = self.closest_frames(time);
        let mut bitmap = vec![false; self.width * self.height];
        if before == after {
            let pixels_len = pixels.len();
            for i in self.keyframes[before]
                .pixels
                .iter()
                .filter(|&i| i < pixels_len)
            {
                pixels[i] = self.color;
            }
            if self.fill_color[3] > 0 {
                for i in self.keyframes[before].fill_pixels.iter() {
                    pixels[i] = self.fill_color;
                }
            }
        } else {
            if self.tweens.get(&(before, after)).is_none() {
                self.tweens.insert(
                    (before, after),
                    Tween::new(
                        self.width,
                        self.keyframes[before].pixels.clone(),
                        self.keyframes[after].pixels.clone(),
                    ),
                );
            }
            let tween = self.tweens.get_mut(&(before, after)).unwrap();
            let fraction = (time - self.keyframes[before].time)
                / (self.keyframes[after].time - self.keyframes[before].time);
            tween.draw(fraction, &mut bitmap);
            for (_, out) in bitmap.iter().zip(pixels.iter_mut()).filter(|(b, _)| **b) {
                *out = self.color;
            }
            if self.fill_color[3] > 0 {
                let mut outside = bitmap.clone();
                let mut todo = vec![0];
                let w = self.width;
                outside[0] = true;
                while let Some(i) = todo.pop() {
                    if i > 0 && !outside[i - 1] {
                        outside[i - 1] = true;
                        todo.push(i - 1);
                    }
                    if i + 1 < outside.len() && !outside[i + 1] {
                        outside[i + 1] = true;
                        todo.push(i + 1);
                    }
                    if i >= w && !outside[i - w] {
                        outside[i - w] = true;
                        todo.push(i - w);
                    }
                    if i + w < outside.len() && !outside[i + w] {
                        outside[i + w] = true;
                        todo.push(i + w);
                    }
                }
                for (_, out) in outside.iter().zip(pixels.iter_mut()).filter(|(b, _)| !**b) {
                    *out = self.fill_color;
                }
            }
        }
    }

    pub fn draw_gif(&mut self, time: f32, color: u8, pixels: &mut [u8]) {
        let (before, after) = self.closest_frames(time);
        let mut bitmap = vec![false; self.width * self.height];
        if before == after {
            let pixels_len = pixels.len();
            for i in self.keyframes[before]
                .pixels
                .iter()
                .filter(|&i| i < pixels_len)
            {
                pixels[i] = color;
            }
            if self.fill_color[3] > 0 {
                for i in self.keyframes[before].fill_pixels.iter() {
                    pixels[i] = color + 1;
                }
            }
        } else {
            if self.tweens.get(&(before, after)).is_none() {
                self.tweens.insert(
                    (before, after),
                    Tween::new(
                        self.width,
                        self.keyframes[before].pixels.clone(),
                        self.keyframes[after].pixels.clone(),
                    ),
                );
            }
            let tween = self.tweens.get_mut(&(before, after)).unwrap();
            let fraction = (time - self.keyframes[before].time)
                / (self.keyframes[after].time - self.keyframes[before].time);
            tween.draw(fraction, &mut bitmap);
            for (_, out) in bitmap.iter().zip(pixels.iter_mut()).filter(|(b, _)| **b) {
                *out = color;
            }
            if self.fill_color[3] > 0 {
                let mut outside = bitmap.clone();
                let mut todo = vec![0];
                let w = self.width;
                outside[0] = true;
                while let Some(i) = todo.pop() {
                    if i > 0 && !outside[i - 1] {
                        outside[i - 1] = true;
                        todo.push(i - 1);
                    }
                    if i + 1 < outside.len() && !outside[i + 1] {
                        outside[i + 1] = true;
                        todo.push(i + 1);
                    }
                    if i >= w && !outside[i - w] {
                        outside[i - w] = true;
                        todo.push(i - w);
                    }
                    if i + w < outside.len() && !outside[i + w] {
                        outside[i + w] = true;
                        todo.push(i + w);
                    }
                }
                for (_, out) in outside.iter().zip(pixels.iter_mut()).filter(|(b, _)| !**b) {
                    *out = color + 1;
                }
            }
        }
    }

    pub fn closest_time(&self, time: f32) -> f32 {
        self.keyframes[self.closest_frame(time)].time
    }
    fn closest_frames(&self, time: f32) -> (usize, usize) {
        let mut above = 0;
        let mut below = 0;

        for (i, f) in self.keyframes.iter().enumerate() {
            if f.time >= time
                && (f.time <= self.keyframes[above].time || self.keyframes[above].time < time)
            {
                above = i;
            }
            if f.time <= time
                && (f.time >= self.keyframes[below].time || self.keyframes[below].time > time)
            {
                below = i;
            }
        }
        if self.keyframes[below].time > time {
            below = above;
        }
        if self.keyframes[above].time < time {
            above = below;
        }
        (below, above)
    }
    fn closest_frame(&self, time: f32) -> usize {
        let mut closest = 2.0;
        for f in self.keyframes.iter() {
            if (f.time - time).abs() < closest {
                closest = (f.time - time).abs();
            }
        }
        if let Some((i, _)) = self
            .keyframes
            .iter()
            .enumerate()
            .filter(|(_, k)| (k.time - time).abs() == closest)
            .next()
        {
            i
        } else {
            0
        }
    }

    pub fn erase_pixels(&mut self, time: f32, pixels: SetUsize) {
        let i = self.closest_frame(time);
        for p in pixels {
            self.keyframes[i].pixels.remove(p);
        }
        self.compute_fill(i);
    }
    pub fn add_pixels(&mut self, time: f32, pixels: SetUsize) {
        let i = self.closest_frame(time);
        self.keyframes[i].pixels.extend(pixels);
        self.compute_fill(i);
    }
    pub fn get_chunk(&self, time: f32, selected: &SetUsize) -> SetUsize {
        let w = self.width as usize;
        let i = self.closest_frame(time);
        let mut pixels = self.keyframes[i].pixels.clone();
        let mut todo: Vec<usize> = selected.iter().filter(|&p| pixels.remove(p)).collect();
        let mut out = SetUsize::new();

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
        out
    }
    pub fn get_filled_chunk(&self, time: f32, selected: &SetUsize) -> SetUsize {
        let w = self.width as usize;
        let i = self.closest_frame(time);
        let mut pixels = self.keyframes[i].pixels.clone();
        pixels.extend(self.keyframes[i].fill_pixels.iter());
        let mut todo: Vec<usize> = selected.iter().collect();
        let mut out = SetUsize::new();
        for p in todo.iter().copied() {
            pixels.remove(p);
        }

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
        out
    }
    pub fn move_pixels(&mut self, time: f32, displacement: Vec2) {
        let i = self.closest_frame(time);
        let didx =
            displacement.x.round() as isize + displacement.y.round() as isize * self.width as isize;
        self.keyframes[i].pixels = self.keyframes[i]
            .pixels
            .iter()
            .map(|idx| idx.wrapping_add(didx as usize))
            .collect();
    }

    pub fn shift_frame(&mut self, old_time: f32, new_time: f32) {
        for (which, k) in self.keyframes.iter_mut().enumerate() {
            if k.time == old_time {
                k.time = new_time;
                self.tweens.retain(|k, _| k.0 != which && k.1 != which);
            }
        }
    }

    pub fn ensure_we_have_frame_at(&mut self, time: f32) {
        let times = self.closest_frames(time);
        if self.keyframes[times.0].time == time || self.keyframes[times.1].time == time {
            // We already have a frame at this time.
            return;
        }
        let mut pixels = SetUsize::new();
        let mut fill_pixels = Pixels::default();
        if times.0 == times.1 {
            // In this case, self.image didn't get updated!
            pixels = self.keyframes[times.0].pixels.clone();
            fill_pixels = self.keyframes[times.0].fill_pixels.clone();
        }
        self.keyframes.push(Bitmap {
            time,
            // texture: Texture2D::from_image(&self.image),
            pixels,
            fill_pixels,
        });
        self.handle_modified_bitmap(time);
    }
}

fn random_color() -> [u8; 4] {
    [rand::random(), rand::random(), rand::random(), 255]
}

pub fn clamp(min: f32, max: f32, val: f32) -> f32 {
    if val < min {
        min
    } else if val > max {
        max
    } else {
        val
    }
}
