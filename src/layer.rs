use std::{
    collections::HashMap,
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
};

use macroquad::{
    prelude::{
        draw_line, draw_rectangle, draw_rectangle_lines, is_mouse_button_down,
        is_mouse_button_pressed, is_mouse_button_released, mouse_position, screen_height,
        screen_width, Color, Image, MouseButton, Texture2D, Vec2, BLACK, GRAY, WHITE,
    },
    texture::{draw_texture_ex, DrawTextureParams},
};
use tinyset::SetUsize;

use crate::tween::Tween;

#[derive(Clone)]
struct Bitmap {
    time: f32,
    pixels: SetUsize,
    fill_pixels: SetUsize,
    texture: Texture2D,
}

pub struct Layer {
    pub color: [u8; 4],
    pub fill_color: [u8; 4],
    image: Image,
    texture: Texture2D,
    keyframes: Vec<Bitmap>,
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
            image: bitmap.clone(),
            texture: Texture2D::from_image(&bitmap),
            keyframes: vec![Bitmap {
                time,
                texture: Texture2D::from_image(&bitmap),
                pixels: SetUsize::new(),
                fill_pixels: SetUsize::new(),
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
    pub fn handle_modified_bitmap(&mut self, time: f32) {
        let which = self.closest_frame(time);
        self.tweens.retain(|k, _| k.0 != which && k.1 != which);
        let k = &mut self.keyframes[which];
        let mut img = Image::gen_image_color(
            self.image.width,
            self.image.height,
            Color::from_rgba(0, 0, 0, 0),
        );
        {
            let img = img.get_image_data_mut();
            let mut outside = vec![false; img.len()];
            let image_len = img.len();
            for i in k.pixels.iter().filter(|&i| i < image_len) {
                img[i] = self.color;
                outside[i] = true;
            }

            let mut todo = vec![0];
            let w = self.image.width();
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
            k.fill_pixels = SetUsize::from_iter(
                outside
                    .into_iter()
                    .enumerate()
                    .filter(|(_, b)| !*b)
                    .map(|x| x.0),
            );
            if self.fill_color[3] > 0 {
                for p in k.fill_pixels.iter() {
                    img[p] = self.fill_color;
                }
            }
        }
        k.texture.update(&img);
    }
    pub fn draw(&mut self, time: f32, pixels: &mut [[u8; 4]]) {
        let (before, after) = self.closest_frames(time);
        let mut bitmap = vec![false; self.image.get_image_data().len()];
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
                        self.image.width(),
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
                let w = self.image.width();
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
    }
    pub fn add_pixels(&mut self, time: f32, pixels: SetUsize) {
        let i = self.closest_frame(time);
        self.keyframes[i].pixels.extend(pixels);
    }
    pub fn move_pixels(&mut self, time: f32, displacement: Vec2) {
        let i = self.closest_frame(time);
        let didx = displacement.x.round() as isize
            + displacement.y.round() as isize * self.image.width() as isize;
        self.keyframes[i].pixels = self.keyframes[i]
            .pixels
            .iter()
            .map(|idx| idx.wrapping_add(didx as usize))
            .collect();
    }

    pub fn frame_selector(&mut self, now: &mut f32) -> bool {
        const TSTART: f32 = 100.0;
        const THEIGHT: f32 = 50.0;
        const FRAME_WIDTH: f32 = 70.0;
        let tstop: f32 = screen_width() - 2.0 * FRAME_WIDTH;
        let t_width = tstop - TSTART;
        draw_line(TSTART, THEIGHT, tstop, THEIGHT, 4.0, GRAY);
        draw_line(
            TSTART + *now * (tstop - TSTART),
            THEIGHT * 0.5,
            TSTART + *now * (tstop - TSTART),
            THEIGHT * 1.5,
            10.0,
            WHITE,
        );

        for frame in self.keyframes.iter() {
            let x = TSTART + frame.time * t_width;
            draw_rectangle(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, BLACK);
            draw_texture_ex(
                frame.texture,
                x,
                THEIGHT * 0.5,
                self.get_color(),
                DrawTextureParams {
                    dest_size: Some(Vec2::new(FRAME_WIDTH, THEIGHT)),
                    ..Default::default()
                },
            );
            let (thickness, color) = if frame.time == *now {
                (4.0, WHITE)
            } else {
                (2.0, GRAY)
            };
            draw_rectangle_lines(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, thickness, color);
        }
        let pos = mouse_position();
        static AM_DRAGGING: AtomicBool = AtomicBool::new(false);
        static KEYFRAME: AtomicUsize = AtomicUsize::new(0);
        let am_dragging = AM_DRAGGING.load(Ordering::Relaxed);
        let drag_frame = KEYFRAME.load(Ordering::Relaxed);
        if am_dragging
            || (pos.1 <= 1.5 * THEIGHT && pos.0 >= TSTART + FRAME_WIDTH * 0.5 && pos.0 < tstop)
        {
            let time = clamp(
                0.0,
                1.0,
                (pos.0 - FRAME_WIDTH * 0.5 - TSTART) / (tstop - TSTART),
            );
            let mouse_released = is_mouse_button_released(MouseButton::Left);
            let mouse_pressed = is_mouse_button_pressed(MouseButton::Left);
            let mouse_down = is_mouse_button_down(MouseButton::Left);
            if let Some((i, frame)) = self
                .keyframes
                .iter()
                .enumerate()
                .filter(|(_, f)| {
                    let x = TSTART + f.time * (tstop - TSTART);
                    (pos.0 >= x - FRAME_WIDTH * 0.5) && (pos.0 <= x + FRAME_WIDTH * 1.5)
                })
                .next()
            {
                let x = TSTART + frame.time * t_width;
                draw_rectangle_lines(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, 4.0, WHITE);
                if am_dragging && drag_frame == i {
                    draw_rectangle_lines(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, 6.0, WHITE);
                    if mouse_released {
                        *now = frame.time;
                        AM_DRAGGING.store(false, Ordering::Relaxed);
                    }
                } else if mouse_pressed {
                    AM_DRAGGING.store(true, Ordering::Relaxed);
                    KEYFRAME.store(i, Ordering::Relaxed);
                }
            } else if am_dragging && self.keyframes.len() > drag_frame {
                let x = clamp(TSTART, tstop, pos.0) - FRAME_WIDTH * 0.5;
                if mouse_down {
                    draw_rectangle(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, BLACK);
                    draw_texture_ex(
                        self.keyframes[drag_frame].texture,
                        x,
                        THEIGHT * 0.5,
                        self.get_color(),
                        DrawTextureParams {
                            dest_size: Some(Vec2::new(FRAME_WIDTH, THEIGHT)),
                            ..Default::default()
                        },
                    );
                    draw_rectangle_lines(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, 2.0, GRAY);
                } else {
                    self.keyframes[KEYFRAME.load(Ordering::Relaxed)].time = time;
                    *now = time;
                    AM_DRAGGING.store(false, Ordering::Relaxed);
                }
            } else {
                draw_rectangle(
                    pos.0 - FRAME_WIDTH * 0.5,
                    THEIGHT * 0.5,
                    FRAME_WIDTH,
                    THEIGHT,
                    BLACK,
                );
                draw_rectangle_lines(
                    pos.0 - FRAME_WIDTH * 0.5,
                    THEIGHT * 0.5,
                    FRAME_WIDTH,
                    THEIGHT,
                    2.0,
                    GRAY,
                );
                draw_line(pos.0 - 10.0, THEIGHT, pos.0 + 10.0, THEIGHT, 4.0, GRAY);
                draw_line(pos.0, THEIGHT - 10.0, pos.0, THEIGHT + 10.0, 4.0, GRAY);
                if mouse_released {
                    let mut img = self.image.clone();
                    self.draw(time, img.get_image_data_mut());
                    self.texture.update(&img);

                    let mut pixels = SetUsize::new();
                    let mut fill_pixels = SetUsize::new();
                    let times = self.closest_frames(time);
                    if times.0 == times.1 {
                        // In this case, self.image didn't get updated!
                        pixels = self.keyframes[times.0].pixels.clone();
                        fill_pixels = self.keyframes[times.0].fill_pixels.clone();
                    }
                    self.keyframes.push(Bitmap {
                        time,
                        texture: Texture2D::from_image(&self.image),
                        pixels,
                        fill_pixels,
                    });
                    self.handle_modified_bitmap(time);
                    *now = time;
                }
            }
        }
        pos.1 < 2.0 * THEIGHT || am_dragging
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
