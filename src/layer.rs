use std::collections::HashMap;

use macroquad::prelude::{screen_height, screen_width, Color, Image, Vec2};
use serde::{Deserialize, Serialize};

use crate::{pixeltree::Pixels, tween::Tween};

#[derive(Clone, Serialize, Deserialize)]
struct Bitmap {
    time: f32,
    pixels: Pixels,
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
                pixels: Pixels::default(),
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
        k.fill_pixels = k.pixels.compute_fill(self.width);
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
                    if let Some(i) = img.get_mut(p) {
                        *i = self.fill_color;
                    }
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
                    if let Some(p) = pixels.get_mut(i) {
                        *p = self.fill_color;
                    }
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

    pub fn draw_points(&mut self, time: f32) -> (Pixels, Pixels) {
        let (before, after) = self.closest_frames(time);
        if before == after {
            let penpoints = self.keyframes[before].pixels.clone();
            let fillpoints = self.keyframes[before].fill_pixels.clone();
            (penpoints, fillpoints)
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
            let points = tween.draw_points(fraction);
            let fill = points.compute_fill(self.width);
            (points, fill)
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
                    if let Some(pix) = pixels.get_mut(i) {
                        *pix = color + 1;
                    }
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
    pub fn closest_points(&self, time: f32) -> &Pixels {
        let i = self.closest_frame(time);
        &self.keyframes[i].pixels
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

    pub fn erase_pixels(&mut self, time: f32, pixels: &Pixels) {
        let i = self.closest_frame(time);
        self.keyframes[i].pixels.remove(pixels);
        self.compute_fill(i);
    }
    pub fn add_pixels(&mut self, time: f32, pixels: &Pixels) {
        let i = self.closest_frame(time);
        self.keyframes[i].pixels.extend(pixels);
        self.compute_fill(i);
    }
    pub fn add_base_rectangle(&mut self, rect_width: usize, rect_height: usize) {
        let mut pixels = Pixels::default();
        let xmin = self.width / 2 - rect_width / 2 - 1;
        let ymin = self.height / 2 - rect_height / 2 - 1;
        let xmax = self.width / 2 + rect_width / 2 + 1;
        let ymax = self.height / 2 + rect_height / 2 + 1;
        for x in xmin..xmax + 1 {
            pixels.insert(x + ymin * self.width);
            pixels.insert(x + ymax * self.width);
        }
        for y in ymin..ymax + 1 {
            pixels.insert(xmin + y * self.width);
            pixels.insert(xmax + y * self.width);
        }
        self.keyframes[0].pixels.extend(&pixels);
        self.compute_fill(0);
    }
    pub fn add_circle(&mut self, diameter: isize) {
        let mut pixels = Pixels::default();
        let centerx = self.width as isize / 2;
        let centery = self.height as isize / 2;
        let xmin = centerx - diameter / 2 - 5;
        let ymin = centery - diameter / 2 - 5;
        let r2min = (diameter / 2).pow(2);
        let r2max = (diameter / 2 + 4).pow(2);
        for x in xmin..(xmin + diameter + 10) {
            for y in ymin..(ymin + diameter + 10) {
                let rad2 = (x - centerx).pow(2) + (y - centery).pow(2);
                if rad2 > r2min && rad2 < r2max {
                    pixels.insert(x as usize + y as usize * self.width);
                }
            }
        }
        self.keyframes[0].pixels.extend(&pixels);
        self.compute_fill(0);
    }
    pub fn get_chunk(&mut self, time: f32, selected: &Pixels) -> Pixels {
        let i = self.closest_frame(time);
        if self.keyframes[i].time != time {
            let (points, fill) = self.draw_points(time);
            if points.borders(selected) || fill.borders(selected) {
                self.ensure_we_have_frame_at(time);
                self.get_chunk(time, selected)
            } else {
                Pixels::default()
            }
        } else {
            self.keyframes[i]
                .pixels
                .contiguous_with(self.width, selected)
        }
    }
    pub fn get_filled_chunk(&self, time: f32, selected: &Pixels) -> Pixels {
        let i = self.closest_frame(time);
        let mut pixels = self.keyframes[i].pixels.clone();
        pixels.extend(&self.keyframes[i].fill_pixels);
        pixels.contiguous_with(self.width, selected)
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
        let (pixels, fill_pixels) = self.draw_points(time);
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
    [
        macroquad::rand::rand() as u8,
        macroquad::rand::rand() as u8,
        macroquad::rand::rand() as u8,
        255,
    ]
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
