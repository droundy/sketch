use macroquad::prelude::{
    draw_circle, draw_circle_lines, draw_line, draw_rectangle, draw_rectangle_lines, draw_texture,
    is_key_pressed, is_mouse_button_down, is_mouse_button_pressed, mouse_position, next_frame,
    screen_height, screen_width, Color, Conf, Image, KeyCode, MouseButton, Texture2D, Vec2, BLACK,
    GRAY, WHITE,
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
}

fn random_color() -> Color {
    Color {
        r: rand::random(),
        g: rand::random(),
        b: rand::random(),
        a: 1.0,
    }
}
