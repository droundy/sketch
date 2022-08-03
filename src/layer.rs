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
}

pub struct Layer {
    pub color: Color,
    keyframes: Vec<Bitmap>,
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
            keyframes: vec![Bitmap {
                time,
                texture: Texture2D::from_image(&bitmap),
                bitmap,
            }],
        }
    }
    pub fn texture(&self, time: f32) -> Texture2D {
        let k = self.closest_frame(time);
        k.texture.update(&k.bitmap);
        k.texture
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
        if pos.1 <= 1.5 * THEIGHT && pos.0 >= TSTART + FRAME_WIDTH * 0.5 {
            let time = (pos.0 - FRAME_WIDTH * 0.5 - TSTART) / (tstop - TSTART);
            static AM_DRAGGING: AtomicBool = AtomicBool::new(false);
            static KEYFRAME: AtomicUsize = AtomicUsize::new(0);
            let am_dragging = AM_DRAGGING.load(Ordering::Relaxed);
            let drag_frame = KEYFRAME.load(Ordering::Relaxed);
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
                if mouse_down {
                    draw_rectangle(
                        pos.0 - FRAME_WIDTH * 0.5,
                        THEIGHT * 0.5,
                        FRAME_WIDTH,
                        THEIGHT,
                        BLACK,
                    );
                    draw_texture_ex(
                        self.keyframes[drag_frame].texture,
                        pos.0 - FRAME_WIDTH * 0.5,
                        THEIGHT * 0.5,
                        self.color,
                        DrawTextureParams {
                            dest_size: Some(Vec2::new(FRAME_WIDTH, THEIGHT)),
                            ..Default::default()
                        },
                    );
                    draw_rectangle_lines(
                        pos.0 - FRAME_WIDTH * 0.5,
                        THEIGHT * 0.5,
                        FRAME_WIDTH,
                        THEIGHT,
                        2.0,
                        GRAY,
                    );
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
                    let bytes = self
                        .closest_frame(time)
                        .bitmap
                        .get_image_data()
                        .iter()
                        .flat_map(|d| d.into_iter())
                        .copied()
                        .collect();
                    let bitmap = Image {
                        bytes,
                        width: self.keyframes[0].bitmap.width,
                        height: self.keyframes[0].bitmap.height,
                    };
                    self.keyframes.push(Bitmap {
                        time,
                        texture: Texture2D::from_image(&bitmap),
                        bitmap,
                    });
                }
            }
        }
        pos.1 < 2.0 * THEIGHT
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
