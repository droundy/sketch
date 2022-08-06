use std::{
    collections::{HashMap, HashSet},
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

#[derive(Clone)]
struct Bitmap {
    time: f32,
    bitmap: Image,
    texture: Texture2D,
    center: Vec2,
}

pub struct Layer {
    pub color: Color,
    image: Image,
    texture: Texture2D,
    keyframes: Vec<Bitmap>,
    connections: HashMap<(usize, usize), Vec<(usize, usize)>>,
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
            image: bitmap.clone(),
            texture: Texture2D::from_image(&bitmap),
            keyframes: vec![Bitmap {
                time,
                center: Vec2::ZERO,
                texture: Texture2D::from_image(&bitmap),
                bitmap,
            }],
            connections: HashMap::new(),
        }
    }
    pub fn handle_modified_bitmap(&mut self, time: f32) {
        let which = self.closest_frame(time);
        self.connections.retain(|k, _| k.0 != which && k.1 != which);
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
        k.texture.update(&k.bitmap);
    }
    pub fn texture(&mut self, time: f32) -> Texture2D {
        let (before, after) = self.closest_frames(time);
        if before == after {
            self.keyframes[before].texture
        } else {
            let w = self.image.width();
            if self.connections.get(&(before, after)).is_none() {
                let mut connections = Vec::new();
                let mut b = HashSet::new();
                let mut a = HashSet::new();
                for (i, _) in self.keyframes[before]
                    .bitmap
                    .get_image_data()
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v[3] == 255)
                {
                    b.insert((i % w, i / w));
                }
                for (i, _) in self.keyframes[after]
                    .bitmap
                    .get_image_data()
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v[3] == 255)
                {
                    a.insert((i % w, i / w));
                }
                let mut b_all_done = false;
                let mut b_unchosen: Vec<_> = b.iter().copied().collect();
                for &(ax, ay) in a.iter() {
                    if let Some((bx, by)) = b_unchosen.pop() {
                        connections.push((bx + by * w, ax + ay * w));
                    } else {
                        b_all_done = true;
                        b_unchosen = b.iter().copied().collect();
                    }
                }
                if !b_all_done {
                    for &(ax, ay) in a.iter() {
                        if let Some((bx, by)) = b_unchosen.pop() {
                            connections.push((bx + by * w, ax + ay * w));
                        } else {
                            break;
                        }
                    }
                }
                self.connections.insert((before, after), connections);
            }
            if let Some(connections) = self.connections.get(&(before, after)) {
                let w_before = (self.keyframes[after].time - time)
                    / (self.keyframes[after].time - self.keyframes[before].time);
                let w_after = (time - self.keyframes[before].time)
                    / (self.keyframes[after].time - self.keyframes[before].time);
                for v in self.image.get_image_data_mut().iter_mut() {
                    *v = [255, 255, 255, 0];
                }
                for &(b, a) in connections.iter() {
                    let bx = b % self.image.width();
                    let by = b / self.image.width();
                    let ax = a % self.image.width();
                    let ay = a / self.image.width();
                    let x = w_before * bx as f32 + w_after * ax as f32;
                    let y = w_before * by as f32 + w_after * ay as f32;
                    let idx = x.round() as usize + (y.round() as usize) * self.image.width();
                    self.image.get_image_data_mut()[idx][3] = 255;
                }
                self.texture.update(&self.image);
                self.texture
            } else {
                unreachable!()
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
                self.color,
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
                        self.color,
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
                    self.texture(time);
                    let bitmap = self.image.clone();
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

fn random_color() -> Color {
    Color {
        r: rand::random(),
        g: rand::random(),
        b: rand::random(),
        a: 1.0,
    }
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
