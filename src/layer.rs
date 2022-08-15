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

use crate::tween::Tween;

#[derive(Clone)]
struct Bitmap {
    time: f32,
    bitmap: Image,
    texture: Texture2D,
    center: Vec2,
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
                center: Vec2::ZERO,
                texture: Texture2D::from_image(&bitmap),
                bitmap,
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
        let mut center = Vec2::ZERO;
        let mut num = 0;
        for (idx, b) in k.bitmap.get_image_data().iter().enumerate() {
            if b[3] > 0 {
                let i = idx % k.bitmap.width();
                let j = idx / k.bitmap.width();
                center += Vec2::new(i as f32, j as f32);
                num += 1;
            }
        }
        k.center = center / num as f32;
        let mut img = k.bitmap.clone();
        for p in img.get_image_data_mut() {
            if p[3] > 0 {
                *p = self.color;
            }
        }
        k.texture.update(&img);
    }
    pub fn draw(&mut self, time: f32, pixels: &mut [[u8; 4]]) {
        let (before, after) = self.closest_frames(time);
        let mut bitmap = vec![false; self.keyframes[before].bitmap.get_image_data().len()];
        if before == after {
            for (p, out) in self.keyframes[before]
                .bitmap
                .get_image_data()
                .iter()
                .zip(bitmap.iter_mut())
            {
                if p[3] != 0 {
                    *out = true;
                }
            }
        } else {
            if self.tweens.get(&(before, after)).is_none() {
                let before_pixels = self.keyframes[before]
                    .bitmap
                    .get_image_data()
                    .iter()
                    .map(|x| x[3] > 0)
                    .collect::<Vec<_>>();
                let after_pixels = self.keyframes[after]
                    .bitmap
                    .get_image_data()
                    .iter()
                    .map(|x| x[3] > 0)
                    .collect::<Vec<_>>();

                self.tweens.insert(
                    (before, after),
                    Tween::new(self.image.width(), before_pixels, after_pixels),
                );
            }
            let tween = self.tweens.get_mut(&(before, after)).unwrap();
            let fraction = (time - self.keyframes[before].time)
                / (self.keyframes[after].time - self.keyframes[before].time);
            tween.draw_bool(fraction, &mut bitmap);
        }
        for (_, out) in bitmap.iter().zip(pixels).filter(|(b,_)| **b) {
            *out = self.color;
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

    pub fn get_frame_data_mut(&mut self, time: f32) -> &mut [[u8; 4]] {
        let i = self.closest_frame(time);
        self.keyframes[i].bitmap.get_image_data_mut()
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

                    let mut bitmap = self.image.clone();
                    let times = self.closest_frames(time);
                    if times.0 == times.1 {
                        // In this case, self.image didn't get updated!
                        bitmap = self.keyframes[times.0].bitmap.clone();
                    }
                    self.keyframes.push(Bitmap {
                        time,
                        center: Vec2::ZERO,
                        texture: Texture2D::from_image(&bitmap),
                        bitmap,
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
