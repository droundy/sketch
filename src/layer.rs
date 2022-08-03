use macroquad::{
    prelude::{
        draw_circle, draw_circle_lines, draw_line, draw_rectangle, draw_rectangle_lines,
        draw_texture, is_key_pressed, is_mouse_button_down, is_mouse_button_pressed,
        mouse_position, next_frame, screen_height, screen_width, Color, Conf, Image, KeyCode,
        MouseButton, Texture2D, Vec2, BLACK, GRAY, WHITE,
    },
    texture::{draw_texture_ex, DrawTextureParams},
};

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
        self.keyframes[0].texture.update(&self.keyframes[0].bitmap);
        self.keyframes[0].texture
    }

    pub fn get_frame_data_mut(&mut self, time: f32) -> &mut [[u8; 4]] {
        self.keyframes[0].bitmap.get_image_data_mut()
    }

    pub fn frame_selector(&mut self) -> bool {
        const TSTART: f32 = 100.0;
        const THEIGHT: f32 = 50.0;
        const FRAME_WIDTH: f32 = 70.0;
        let tstop: f32 = screen_width() - 2.0 * FRAME_WIDTH;
        let t_width = tstop - TSTART;
        draw_line(TSTART, THEIGHT, tstop, THEIGHT, 4.0, GRAY);

        for frame in self.keyframes.iter() {
            let time = frame.time;
            let x = TSTART + time * t_width;
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
            draw_rectangle_lines(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, 2.0, WHITE);
        }
        let pos = mouse_position();
        if pos.1 < 2.0 * THEIGHT {
            let time = (pos.0 - TSTART) / (tstop - TSTART);
            if let Some(frame) = self
                .keyframes
                .iter()
                .filter(|f| {
                    let x = TSTART + f.time * (tstop - TSTART);
                    (pos.0 >= x - FRAME_WIDTH*0.5) && (pos.0 <= x + FRAME_WIDTH)
                })
                .next()
            {
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
