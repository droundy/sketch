use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

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
}

fn shift_img(width: usize, output: &mut [[u8; 4]], input: &[[u8; 4]], shift: Vec2, w: f32) {
    let offset = shift.x.round() + shift.y.round() * width as f32;
    if offset < 0.0 {
        let offset = (-offset) as usize;
        for (i, o) in input[offset..]
            .iter()
            .zip(output[..input.len() - offset].iter_mut())
        {
            o[3] += (w * i[3] as f32).round() as u8;
        }
    } else {
        let offset = offset as usize;
        for (i, o) in input[..input.len() - offset]
            .iter()
            .zip(output[offset..].iter_mut())
        {
            o[3] += (w * i[3] as f32).round() as u8;
        }
    }
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
        }
    }
    pub fn handle_modified_bitmap(&mut self, time: f32) {
        let k = self.closest_frame_mut(time);
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
        let width = before.bitmap.width();
        let mut w_before = 1.0;
        let mut w_after = 0.0;
        if after.time > before.time {
            w_before = (after.time - time) / (after.time - before.time);
            w_after = (time - before.time) / (after.time - before.time);
        }
        let mut out = vec![[255, 255, 255, 0]; before.bitmap.get_image_data().len()];
        let center = w_before * before.center + w_after * after.center;
        shift_img(
            width,
            &mut out,
            before.bitmap.get_image_data(),
            center - before.center,
            w_before,
        );
        shift_img(
            width,
            &mut out,
            after.bitmap.get_image_data(),
            center - after.center,
            w_after,
        );
        for (o, i) in self
            .image
            .get_image_data_mut()
            .iter_mut()
            .zip(out.into_iter())
        {
            *o = i;
        }
        self.texture.update(&self.image);
        self.texture
        // let k = self.closest_frame(time);
        // k.texture
    }
    pub fn closest_time(&self, time: f32) -> f32 {
        self.closest_frame(time).time
    }
    fn closest_frames(&self, time: f32) -> (&Bitmap, &Bitmap) {
        let mut above = &self.keyframes[0];
        let mut below = &self.keyframes[0];

        for f in self.keyframes.iter() {
            if f.time >= time && (f.time < above.time || above.time < time) {
                above = f;
            }
            if f.time <= time && (f.time > below.time || below.time > time) {
                below = f;
            }
        }
        (below, above)
    }
    fn closest_frame(&self, time: f32) -> &Bitmap {
        let mut closest = 2.0;
        for f in self.keyframes.iter() {
            if (f.time - time).abs() < closest {
                closest = (f.time - time).abs();
            }
        }
        if let Some(k) = self
            .keyframes
            .iter()
            .filter(|k| (k.time - time).abs() == closest)
            .next()
        {
            k
        } else {
            &self.keyframes[0]
        }
    }
    fn closest_frame_mut(&mut self, time: f32) -> &mut Bitmap {
        let mut closest = 2.0;
        for f in self.keyframes.iter() {
            if (f.time - time).abs() < closest {
                closest = (f.time - time).abs();
            }
        }
        self.keyframes
            .iter_mut()
            .filter(|k| (k.time - time).abs() == closest)
            .next()
            .unwrap()
    }

    pub fn get_frame_data_mut(&mut self, time: f32) -> &mut [[u8; 4]] {
        self.closest_frame_mut(time).bitmap.get_image_data_mut()
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
