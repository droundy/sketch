use std::borrow::Cow;
use std::path::Path;
use std::time::Instant;
use std::{
    f32::consts::PI,
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
};

use macroquad::prelude::{
    draw_circle, draw_circle_lines, draw_line, draw_rectangle, draw_rectangle_lines, draw_texture,
    is_key_pressed, is_mouse_button_down, is_mouse_button_pressed, is_mouse_button_released,
    mouse_position, next_frame, screen_height, screen_width, Color, Conf, KeyCode, MouseButton,
    Vec2, BLACK, DARKGRAY, GRAY, WHITE,
};
use macroquad::shapes::{draw_poly, draw_triangle, draw_triangle_lines};
use macroquad::texture::{draw_texture_ex, DrawTextureParams, Image, Texture2D};
use macroquad::ui::root_ui;

mod layer;
mod pixeltree;
mod tween;
use layer::{clamp, Layer};
use pixeltree::Pixels;
use serde::{Deserialize, Serialize};

fn conf() -> Conf {
    Conf {
        window_title: String::from("Fun draw"),
        window_width: 1920,
        window_height: 1080,
        fullscreen: true,
        ..Default::default()
    }
}

fn color_selector_color(fx: f32, fy: f32) -> Option<[u8; 4]> {
    let fx = fx * 2.0 - 1.0;
    let fy = fy * 2.0 - 1.0;
    let angle = fy.atan2(fx) + PI;
    let radius = (fx * fx + fy * fy).sqrt();
    const PI6: f32 = PI / 3.0;
    let c1 = (((angle + PI6) % (2.0 * PI / 3.0)) / (PI / 3.0) - 1.0).abs();
    let c1 = 1.0 - (1.0 - c1).abs().powf(1.2);
    assert!(c1.abs() <= 1.0);
    let rgb = if angle < PI6 {
        [1.0, c1, 0.0]
    } else if angle < 2.0 * PI6 {
        [c1, 1.0, 0.0]
    } else if angle < 3.0 * PI6 {
        [0.0, 1.0, c1]
    } else if angle < 4.0 * PI6 {
        [0.0, c1, 1.0]
    } else if angle < 5.0 * PI6 {
        [c1, 0.0, 1.0]
    } else {
        [1.0, 0.0, c1]
    };
    const RSATURATED: f32 = 0.6;
    if radius < RSATURATED {
        let x = 1.0 - radius / RSATURATED;
        let extra = x * (1.0 - (x - 1.0).powi(4));
        let rgb = rgb.map(|c| ((c * (1.0 - extra) + extra) * 255.0) as u8);
        Some([rgb[0], rgb[1], rgb[2], 255])
    } else if radius <= 1.0 {
        let x = 1.0 - (radius - RSATURATED) / (1.0 - RSATURATED);
        let r = x * (1.0 - (1.0 - x).powi(2));
        let rgb = rgb.map(|c| (r * c * 255.0) as u8);
        Some([rgb[0], rgb[1], rgb[2], 255])
    } else if fx > 0.0 && fy > 0.0 {
        Some([0, 0, 0, 0])
    } else {
        None
    }
}

impl Drawing {
    fn open(path: &str, images: &mut Vec<Image>, textures: &mut Vec<Texture2D>) -> Option<Self> {
        let bytes = std::fs::read(path).ok()?;
        let mut drawing: Drawing = if path.ends_with(".json") {
            serde_json::from_slice(&bytes).ok()?
        } else {
            bincode::deserialize(&bytes).ok()?
        };
        for f in drawing.keyframes.iter() {
            for l in drawing.layers.iter_mut() {
                l.handle_modified_bitmap(f.time);
            }
        }
        drawing.handle_modified_bitmap(images, textures);
        Some(drawing)
    }
    fn save(&self, path: &str) -> Result<(), std::io::Error> {
        std::fs::write("drawing.json", serde_json::to_string(self).unwrap())?;
        if path.ends_with(".json") {
            std::fs::write(path, serde_json::to_string(self).unwrap())
        } else {
            std::fs::write(path, bincode::serialize(self).unwrap())
        }
    }
    fn pen_drew(&mut self, pixels: Pixels) {
        if self.keyframes.iter().any(|k| k.time == self.time) {
            self.layers[self.current].ensure_we_have_frame_at(self.time);
        }
        if self.tool == Tool::Eraser {
            self.layers[self.current].erase_pixels(self.time, &pixels);
        } else {
            self.layers[self.current].add_pixels(self.time, &pixels);
        }
    }
    fn move_pixels(&mut self, displacement: Vec2) {
        self.layers[self.current].ensure_we_have_frame_at(self.time);
        self.layers[self.current].move_pixels(self.time, displacement);
    }
    fn move_chunk(
        &mut self,
        old_position: Vec2,
        new_position: Vec2,
        moving_chunk: &mut MovingChunk,
    ) {
        let w = self.width as usize;
        if moving_chunk.the_mask.is_empty() {
            self.layers[self.current].ensure_we_have_frame_at(self.time);
            let idx = old_position.x.round() as usize + old_position.y.round() as usize * w;
            println!("Creating moving chunk from {idx}");
            *moving_chunk = MovingChunk::from_mask(
                self.layers[self.current]
                    .get_filled_chunk(self.time, &Pixels::from_iter([idx]))
                    .iter()
                    .collect(),
            );
            for i in self.current..self.layers.len() {
                moving_chunk
                    .insert_chunk(i, self.layers[i].get_chunk(self.time, moving_chunk.mask()));
            }
            println!("Done setting up moving_chunk.");
        }
        let offset = new_position - old_position;
        let offset = offset.x.round() as isize + offset.y.round() as isize * self.width as isize;

        // First we erase the chunk, and check whether the moved chunk is going to collide
        // with another chunk.
        let mut shifted = moving_chunk.clone();
        shifted.shift_by(offset);
        let mut need_to_not_move = false;
        for which in self.current..self.layers.len() {
            let l = &mut self.layers[which];
            if !moving_chunk.chunks(which).is_empty() {
                l.erase_pixels(self.time, moving_chunk.chunks(which));
            }
            let pixels = l.closest_points(self.time);
            if shifted.mask().borders(pixels) || shifted.chunks(which).borders(pixels) {
                need_to_not_move = true;
            }
        }
        if !need_to_not_move {
            // We can offset the chunk!
            *moving_chunk = shifted;
        }
        // Now draw the chunk in its new (or maybe old?) position.
        for which in self.current..self.layers.len() {
            let l = &mut self.layers[which];
            l.add_pixels(self.time, moving_chunk.chunks(which));
        }
    }
    pub fn frame_selector(
        &mut self,
        images: &mut Vec<Image>,
        textures: &mut Vec<Texture2D>,
    ) -> bool {
        const TSTART: f32 = 100.0;
        const THEIGHT: f32 = 50.0;
        const FRAME_WIDTH: f32 = 70.0;
        let tstop: f32 = screen_width() - 2.0 * FRAME_WIDTH;
        let t_width = tstop - TSTART;
        draw_line(TSTART, THEIGHT, tstop, THEIGHT, 4.0, GRAY);
        draw_line(
            TSTART + self.time * (tstop - TSTART),
            THEIGHT * 0.5,
            TSTART + self.time * (tstop - TSTART),
            THEIGHT * 1.5,
            10.0,
            WHITE,
        );

        for (whichframe, frame) in self.keyframes.iter().enumerate() {
            let x = TSTART + frame.time * t_width;
            draw_rectangle(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, BLACK);
            if whichframe >= images.len() {
                images.push(Image::gen_image_color(
                    self.width as u16,
                    self.height as u16,
                    BLACK,
                ))
            }
            if whichframe >= textures.len() {
                textures.push(Texture2D::from_image(&images[whichframe]));
            }
            draw_texture_ex(
                textures[whichframe],
                x,
                THEIGHT * 0.5,
                WHITE,
                DrawTextureParams {
                    dest_size: Some(Vec2::new(FRAME_WIDTH, THEIGHT)),
                    ..Default::default()
                },
            );
            let (thickness, color) = if frame.time == self.time {
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
                        self.time = frame.time;
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
                        textures[drag_frame],
                        x,
                        THEIGHT * 0.5,
                        WHITE,
                        DrawTextureParams {
                            dest_size: Some(Vec2::new(FRAME_WIDTH, THEIGHT)),
                            ..Default::default()
                        },
                    );
                    draw_rectangle_lines(x, THEIGHT * 0.5, FRAME_WIDTH, THEIGHT, 2.0, GRAY);
                    self.time = time;
                } else {
                    let which_frame = KEYFRAME.load(Ordering::Relaxed);
                    let old_time = self.keyframes[which_frame].time;
                    for l in self.layers.iter_mut() {
                        l.shift_frame(old_time, time);
                    }
                    self.keyframes[which_frame].time = time;
                    self.time = time;
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
                    self.keyframes.push(Keyframe { time });
                    self.time = time;
                    self.handle_modified_bitmap(images, textures);
                }
            }
        }
        pos.1 < 2.0 * THEIGHT || am_dragging
    }

    fn handle_modified_bitmap(&mut self, images: &mut Vec<Image>, textures: &mut Vec<Texture2D>) {
        self.layers[self.current].handle_modified_bitmap(self.time);
        for (i, f) in self.keyframes.iter_mut().enumerate() {
            while images.len() <= i {
                images.push(Image::gen_image_color(
                    self.width,
                    self.height,
                    Color::from_rgba(20, 20, 20, 255),
                ));
            }
            images[i] =
                Image::gen_image_color(self.width, self.height, Color::from_rgba(20, 20, 20, 255));
            for l in self.layers.iter_mut() {
                l.draw(f.time, images[i].get_image_data_mut());
            }
            while textures.len() <= i {
                textures.push(Texture2D::from_image(&images[i]));
            }
            textures[i].update(&images[i]);
        }
    }
    fn animation_button(&mut self) -> bool {
        const RADIUS: f32 = 16.0;
        let center = Vec2::new(screen_width() - 80.0, 50.0);
        if !self.am_animating {
            draw_poly(center.x - 3.0, center.y, 3, RADIUS, 0.0, WHITE);
        } else {
            draw_rectangle(center.x - 10.0, center.y - 10.0, 20.0, 20.0, WHITE);
        }
        draw_rectangle_lines(center.x - 30.0, center.y - 25.0, 60.0, 50.0, 4.0, WHITE);
        let (x, y) = mouse_position();
        if is_mouse_button_pressed(MouseButton::Left)
            && (x - center.x).abs() < 30.0
            && (y - center.y).abs() < 25.0
        {
            self.am_animating = !self.am_animating;
            if !self.am_animating && !is_mouse_button_released(MouseButton::Left) {
                self.time = self.layers[self.current].closest_time(self.time);
            }
            true
        } else {
            false
        }
    }
    fn color_selector(&mut self) -> bool {
        let swidth = screen_width();
        let sheight = screen_height();

        let w = 200.0;
        let h = 200.0;
        let dx = w * 0.01;
        let dy = h * 0.01;
        for i in 0..100 {
            for j in 0..100 {
                let fx = i as f32 * 0.01;
                let fy = j as f32 * 0.01;
                if let Some(c) = color_selector_color(fx, fy) {
                    let x = fx * w + swidth - w;
                    let y = fy * h + sheight - h;
                    draw_rectangle(x, y, dx, dy, Color::from_rgba(c[0], c[1], c[2], c[3]));
                }
            }
        }

        static AM_DRAGGING: AtomicBool = AtomicBool::new(false);
        let pos = mouse_position();
        let fx = (pos.0 + w - swidth) / w;
        let fy = (pos.1 + h - sheight) / h;
        if let Some(c) = color_selector_color(fx, fy) {
            if is_mouse_button_pressed(MouseButton::Left) {
                AM_DRAGGING.store(true, std::sync::atomic::Ordering::Relaxed);
            } else if !is_mouse_button_down(MouseButton::Left) {
                AM_DRAGGING.store(false, std::sync::atomic::Ordering::Relaxed);
            }
            if AM_DRAGGING.load(std::sync::atomic::Ordering::Relaxed) {
                if self.am_selecting_fill {
                    self.layers[self.current].fill_color = c;
                } else {
                    self.layers[self.current].color = c;
                }
            }
            true
        } else {
            AM_DRAGGING.load(std::sync::atomic::Ordering::Relaxed)
        }
    }
    fn layer_selector(&mut self) -> bool {
        const WIDTH: f32 = 80.0;
        const HEIGHT: f32 = 80.0;
        let bottom_layer_index = (screen_height() / HEIGHT) as usize - 1;
        let bottom_layer_y = bottom_layer_index as f32 * HEIGHT;
        fn outline(y: f32, selected: bool) {
            draw_rectangle_lines(0.0, y, WIDTH, HEIGHT, 5.0, BLACK);
            draw_rectangle_lines(
                2.0,
                y + 2.0,
                WIDTH - 4.0,
                HEIGHT - 4.0,
                5.0,
                if selected { WHITE } else { GRAY },
            );
        }
        for (i, l) in self.layers.iter().enumerate() {
            let y = bottom_layer_y - i as f32 * HEIGHT;
            draw_rectangle(0.0, y, WIDTH, HEIGHT, l.get_color());
            draw_rectangle(
                0.3 * WIDTH,
                y + 0.3 * HEIGHT,
                WIDTH * 0.4,
                HEIGHT * 0.4,
                if l.fill_color[3] == 0 {
                    DARKGRAY
                } else {
                    l.get_fill_color()
                },
            );
            outline(y, i == self.current && !self.am_selecting_fill);
            if self.am_selecting_fill && i == self.current {
                draw_rectangle_lines(
                    0.3 * WIDTH,
                    y + 0.3 * HEIGHT,
                    WIDTH * 0.4,
                    HEIGHT * 0.4,
                    4.0,
                    WHITE,
                );
            }
        }
        let y = bottom_layer_y - self.layers.len() as f32 * HEIGHT;
        // Draw new layer button
        if y > 3.0 * HEIGHT {
            draw_rectangle(0.0, y, WIDTH, HEIGHT, BLACK);
            outline(y, false);
            draw_line(
                WIDTH * 0.5,
                y + HEIGHT * 0.5 - 10.0,
                WIDTH * 0.5,
                y + HEIGHT * 0.5 + 10.0,
                5.0,
                WHITE,
            );
            draw_line(
                WIDTH * 0.5 - 10.0,
                y + HEIGHT * 0.5,
                WIDTH * 0.5 + 10.0,
                y + HEIGHT * 0.5,
                5.0,
                WHITE,
            );
        }
        // Now draw the tool buttons.
        let y = 0.0;
        outline(y, self.tool == Tool::BigPen);
        self.draw_tool(Tool::BigPen, WIDTH * 0.5, y + HEIGHT * 0.5, true);
        outline(y + HEIGHT, self.tool == Tool::LittlePen);
        self.draw_tool(Tool::LittlePen, WIDTH * 0.5, y + HEIGHT * 1.5, true);
        outline(y + 2.0 * HEIGHT, self.tool == Tool::Eraser);
        self.draw_tool(Tool::Eraser, WIDTH * 0.5, y + HEIGHT * 2.5, true);
        // outline(y + 3.0 * HEIGHT, self.tool == Tool::Move);
        // self.draw_tool(Tool::Move, WIDTH * 0.5, y + HEIGHT * 3.5, true);
        outline(y + 3.0 * HEIGHT, self.tool == Tool::MoveChunk);
        self.draw_tool(Tool::MoveChunk, WIDTH * 0.5, y + HEIGHT * 3.5, true);
        let (x, y) = mouse_position();
        if let Some(original_layer) = self.am_dragging_layer {
            if is_mouse_button_released(MouseButton::Left) {
                self.am_dragging_layer = None;
            } else {
                let idx = (y / HEIGHT) as usize;
                let idx = if idx > bottom_layer_index {
                    0
                } else {
                    bottom_layer_index - idx
                };
                let idx = std::cmp::min(idx, self.layers.len() - 1);
                if idx != original_layer {
                    // Only animate dragging if we have dragged off the original layer.
                    draw_rectangle(
                        0.0,
                        bottom_layer_y - self.current as f32 * HEIGHT,
                        WIDTH,
                        HEIGHT,
                        BLACK,
                    );
                    draw_rectangle(
                        0.0,
                        y - 0.5 * HEIGHT,
                        WIDTH,
                        HEIGHT,
                        self.layers[self.current].get_color(),
                    );
                    draw_rectangle(
                        0.3 * WIDTH,
                        y - 0.2 * HEIGHT,
                        WIDTH * 0.4,
                        HEIGHT * 0.4,
                        if self.layers[self.current].fill_color[3] == 0 {
                            BLACK
                        } else {
                            self.layers[self.current].get_fill_color()
                        },
                    );
                    outline(y - 0.5 * HEIGHT, true);
                }
                if self.current != idx {
                    let t = self.layers.remove(self.current);
                    self.layers.insert(idx, t);
                    self.current = idx;
                }
            }
            return true;
        } else if is_mouse_button_pressed(MouseButton::Left) {
            let y = (y / HEIGHT) as usize;
            if x < WIDTH {
                if y == 0 {
                    self.tool = Tool::BigPen;
                } else if y == 1 {
                    self.tool = Tool::LittlePen;
                } else if y == 2 {
                    self.tool = Tool::Eraser;
                // } else if y == 3 {
                //     self.tool = Tool::Move;
                } else if y == 3 {
                    self.tool = Tool::MoveChunk;
                } else if y == bottom_layer_index - self.layers.len() {
                    self.layers.push(Layer::new(self.time));
                    self.current = self.layers.len() - 1;
                    if self.tool == Tool::Eraser {
                        self.tool = Tool::BigPen;
                    }
                    self.am_selecting_fill = false;
                } else {
                    let new_current = bottom_layer_index - y;
                    if new_current < self.layers.len() {
                        self.am_dragging_layer = Some(new_current);
                        self.current = bottom_layer_index - y;
                        let (x, y) = mouse_position();
                        let y = y % HEIGHT;
                        self.am_selecting_fill = x > 0.3 * WIDTH
                            && x < 0.7 * WIDTH
                            && y > 0.3 * HEIGHT
                            && y < 0.7 * HEIGHT;
                    }
                }
                return true;
            }
        }
        x < WIDTH
    }

    fn show_cursor(&self) {
        let (x, y) = mouse_position();
        self.draw_tool(self.tool, x, y, false);
    }

    fn draw_tool(&self, tool: Tool, x: f32, y: f32, filled: bool) {
        let color = self.layers[self.current].color;
        let contrast = if color[0] as u64 + color[1] as u64 + color[2] as u64 > 300 {
            BLACK
        } else {
            WHITE
        };
        const RADIUS: f32 = Tool::Eraser.radius() * 0.7 + Tool::BigPen.radius() * 0.3;
        let color = Color::from_rgba(color[0], color[1], color[2], color[3]);
        if tool == Tool::Move || tool == Tool::MoveChunk {
            draw_line(x, y - RADIUS, x, y + RADIUS, 4.0, BLACK);
            draw_line(x - RADIUS, y, x + RADIUS, y, 4.0, BLACK);
            draw_line(x, y - RADIUS, x, y + RADIUS, 2.0, WHITE);
            draw_line(x - RADIUS, y, x + RADIUS, y, 2.0, WHITE);
            let c = Vec2::new(x, y);
            for point in [
                Vec2::new(RADIUS, 0.0),
                Vec2::new(-RADIUS, 0.0),
                Vec2::new(0.0, -RADIUS),
                Vec2::new(0.0, RADIUS),
            ] {
                let orth = Vec2::new(point.y, -point.x);
                draw_triangle_lines(
                    c + point,
                    c + point * 0.6 + orth * 0.3,
                    c + point * 0.6 - orth * 0.3,
                    3.0,
                    BLACK,
                );
                draw_triangle(
                    c + point,
                    c + point * 0.6 + orth * 0.3,
                    c + point * 0.6 - orth * 0.3,
                    WHITE,
                );
            }
            if tool == Tool::MoveChunk {
                draw_circle_lines(x, y, RADIUS * 1.1, 2.0, color);
                draw_circle_lines(x, y, RADIUS * 1.1 + 1.0, 1.0, contrast);
            }
        } else if tool == Tool::Eraser {
            draw_circle_lines(x, y, tool.radius(), 1.0, BLACK);
            draw_circle_lines(x, y, tool.radius() + 1.0, 1.0, WHITE);
        } else if tool == Tool::LittlePen {
            draw_circle(x, y, tool.radius(), color);
            draw_circle_lines(x, y, tool.radius() + 1.0, 1.0, contrast);
        } else {
            if filled {
                draw_circle(x, y, tool.radius(), color);
            } else {
                draw_circle_lines(x, y, tool.radius() - 2.0, 3.0, color);
            }
            draw_circle_lines(x, y, tool.radius() + 1.0, 1.0, contrast);
        }
    }
}

#[derive(PartialEq, Eq, Debug, Copy, Clone, Deserialize, Serialize)]
enum Tool {
    Eraser,
    BigPen,
    LittlePen,
    Move,
    MoveChunk,
}

impl Tool {
    const fn radius(&self) -> f32 {
        match self {
            Tool::LittlePen => 2.0,
            Tool::BigPen => 10.0,
            Tool::Eraser => 20.0,
            Tool::Move => 0.0,
            Tool::MoveChunk => 0.0,
        }
    }
}

#[derive(Clone, Deserialize, Serialize)]
struct Drawing {
    am_animating: bool,
    am_dragging_layer: Option<usize>,
    am_selecting_fill: bool,
    current: usize,
    height: u16,
    width: u16,
    time: f32,
    tool: Tool,
    layers: Vec<Layer>,
    keyframes: Vec<Keyframe>,
}

#[derive(Clone, Deserialize, Serialize)]
struct Keyframe {
    time: f32,
}

#[macroquad::main(conf)]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut filename = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "drawing.sketch".to_owned());
    let mut bitmap = Image::gen_image_color(
        screen_width() as u16,
        screen_height() as u16,
        Color {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 0.0,
        },
    );
    let texture = Texture2D::from_image(&bitmap);
    macroquad::input::show_mouse(false);

    let mut old_pos: Option<Vec2> = None;
    let mut frame_images = Vec::new();
    let mut frame_textures = Vec::new();
    let mut drawing =
        Drawing::open(&filename, &mut frame_images, &mut frame_textures).unwrap_or(Drawing {
            am_animating: false,
            am_dragging_layer: None,
            am_selecting_fill: false,
            time: 0.0,
            current: 0,
            tool: Tool::BigPen,
            width: screen_width() as u16,
            height: screen_height() as u16,
            layers: vec![Layer::new(0.0)],
            keyframes: vec![Keyframe { time: 0.0 }],
        });
    let width = drawing.width as usize;
    let height = drawing.height as usize;
    let mut started = Instant::now();
    let mut needs_save = false;
    let mut moving_chunk = MovingChunk::from_mask(Pixels::default());
    loop {
        // clear_background(WHITE);
        if is_key_pressed(KeyCode::Escape)
            || is_key_pressed(KeyCode::Tab)
            || is_key_pressed(KeyCode::Right)
            || is_key_pressed(KeyCode::Left)
            || is_key_pressed(KeyCode::Space)
        {
            if needs_save {
                drawing.save(&filename).ok();
            }
            let mut dir = Path::new(&filename).parent().unwrap_or(Path::new("."));
            std::fs::create_dir_all(dir).expect("Unable to create sketch directory!");
            let gifname = format!("{filename}.gif");
            let mut image = std::fs::File::create(&gifname).unwrap();
            let mut color_map = Vec::with_capacity(drawing.layers.len() * 6);
            color_map.extend([0, 0, 0]);
            for l in drawing.layers.iter() {
                color_map.extend(&l.color[0..3]);
                color_map.extend(&l.fill_color[0..3]);
            }
            let mut encoder =
                gif::Encoder::new(&mut image, width as u16, height as u16, &color_map).unwrap();
            encoder.set_repeat(gif::Repeat::Infinite).unwrap();
            let mut t = 0.0;
            let mut old_pixels = Vec::new();
            let mut delay = 0;
            while t <= 1.0 {
                let mut frame = gif::Frame::default();
                frame.width = width as u16;
                frame.height = height as u16;
                let mut buf = vec![0; width * height];
                for (i, l) in drawing.layers.iter_mut().enumerate() {
                    l.draw_gif(t, 2 * i as u8 + 1, &mut buf);
                }
                frame.buffer = Cow::Borrowed(&buf);
                delay += 2;
                frame.delay = delay;
                // Only write a new frame if something has changed.
                if buf != old_pixels || t + 0.01 > 1.0 {
                    encoder.write_frame(&frame).unwrap();
                    delay = 0;
                }
                print!("\rCreating gif: {t:4.2}");
                t += 0.01;
                old_pixels = buf;
            }
            if is_key_pressed(KeyCode::Escape) {
                return;
            }
            needs_save = false;
            moving_chunk = MovingChunk::from_mask(Pixels::default());
            frame_images = Vec::new();
            for t in frame_textures.drain(..) {
                t.delete();
            }
            let basename = Path::new(&filename)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .to_string();
            println!("looking in directory {dir:?}");
            if dir == Path::new("") {
                dir = Path::new(".");
            }
            let mut files = dir
                .read_dir()
                .unwrap()
                .flat_map(|e| e.ok())
                .flat_map(|e| e.file_name().into_string().ok())
                .filter(|f| f.ends_with(".json"))
                .collect::<Vec<_>>();
            files.sort();
            if files.len() < 2 || is_key_pressed(KeyCode::Space) {
                files = Vec::new();
            }
            if let Ok(i) = files.binary_search(&basename) {
                let order: Vec<_> = if is_key_pressed(KeyCode::Left) {
                    ((0..i).rev()).chain((i + 1..files.len()).rev()).collect()
                } else {
                    (i + 1..files.len()).chain(0..i).collect()
                };
                for i in order {
                    filename = files[i].clone();
                    println!("Opening file {filename}");
                    if let Some(d) =
                        Drawing::open(&filename, &mut frame_images, &mut frame_textures)
                    {
                        if width == d.width as usize && height == d.height as usize {
                            drawing = d;
                            break;
                        } else {
                            println!(
                                "Unable to use file {filename} because it is {} x {} rather than {width} x {height}",
                                d.width, d.height
                            );
                        }
                    }
                }
            } else {
                println!("There is no file {basename} in {dir:?}");
                let mut i = 0;
                let mut f = dir.join(format!("sketch-{i:03x}.json"));
                while f.exists() {
                    i += 1;
                    f = dir.join(format!("sketch-{i:03x}.json"));
                }
                filename = f.to_str().unwrap().to_string();
                drawing = Drawing {
                    am_animating: false,
                    am_dragging_layer: None,
                    am_selecting_fill: false,
                    time: 0.0,
                    current: 0,
                    tool: Tool::BigPen,
                    width: screen_width() as u16,
                    height: screen_height() as u16,
                    layers: vec![Layer::new(0.0)],
                    keyframes: vec![Keyframe { time: 0.0 }],
                };
            }
        }
        if drawing.am_animating {
            drawing.time = (started.elapsed().as_secs_f32() * 0.2) % 1.0;
        }

        for b in bitmap.get_image_data_mut().iter_mut() {
            *b = [20, 20, 20, 255];
        }
        for l in drawing.layers.iter_mut() {
            l.draw(drawing.time, bitmap.get_image_data_mut());
        }
        texture.update(&bitmap);
        draw_texture(texture, 0.0, 0.0, WHITE);

        let animation_button_selected = drawing.animation_button();
        if animation_button_selected {
            started = Instant::now();
        }
        let color_selected = drawing.color_selector();
        if color_selected {
            needs_save = true;
            drawing.handle_modified_bitmap(&mut frame_images, &mut frame_textures);
        }
        let frame_selected = drawing.frame_selector(&mut frame_images, &mut frame_textures);
        let layer_selected = drawing.layer_selector();
        if !root_ui().is_mouse_captured()
            && !color_selected
            && !layer_selected
            && !frame_selected
            && !animation_button_selected
            && !drawing.am_animating
        {
            if is_mouse_button_down(MouseButton::Left) {
                let pos = mouse_position();
                let pos = Vec2::new(pos.0, pos.1);
                let radius = drawing.tool.radius();
                let mut drawn = Pixels::default();
                if let Some(old) = old_pos {
                    if drawing.tool == Tool::Move {
                        drawing.move_pixels(pos - old);
                    } else if drawing.tool == Tool::MoveChunk {
                        drawing.move_chunk(old, pos, &mut moving_chunk);
                    }
                    let parallel = (old - pos).normalize();
                    let orthog = Vec2::new(parallel.y, -parallel.x);
                    let offset = pos.dot(orthog);

                    let mut par_start = pos.dot(parallel);
                    let mut par_stop = old.dot(parallel);
                    if par_stop < par_start {
                        std::mem::swap(&mut par_stop, &mut par_start);
                    }
                    let x_start = std::cmp::max(
                        radius as usize,
                        std::cmp::min(pos.x as usize, old.x as usize),
                    ) - radius as usize;
                    let x_stop = std::cmp::min(
                        width - 1 - (radius as usize),
                        std::cmp::max(pos.x as usize, old.x as usize),
                    ) + radius as usize
                        + 1;
                    let y_start = std::cmp::max(
                        radius as usize,
                        std::cmp::min(pos.y as usize, old.y as usize),
                    ) - radius as usize;
                    let y_stop = std::cmp::min(
                        height - 1 - (radius as usize),
                        std::cmp::max(pos.y as usize, old.y as usize),
                    ) + radius as usize
                        + 1;
                    for x in x_start..x_stop {
                        for y in y_start..y_stop {
                            let here = Vec2::new(x as f32, y as f32);
                            if (here.dot(orthog) - offset).abs() < radius
                                && here.dot(parallel) > par_start
                                && here.dot(parallel) < par_stop
                            {
                                drawn.insert(x + y * width);
                            }
                        }
                    }
                }
                let x_start = std::cmp::max(radius as usize, pos.x as usize) - radius as usize;
                let x_stop = std::cmp::min(width - 1 - (radius as usize), pos.x as usize)
                    + radius as usize
                    + 1;
                let y_start = std::cmp::max(radius as usize, pos.y as usize) - radius as usize;
                let y_stop = std::cmp::min(height - 1 - (radius as usize), pos.y as usize)
                    + radius as usize
                    + 1;
                for x in x_start..x_stop {
                    for y in y_start..y_stop {
                        if (x as f32 - pos.x).powi(2) + (y as f32 - pos.y).powi(2) < radius.powi(2)
                        {
                            drawn.insert(x + y * width);
                        }
                        old_pos = Some(pos);
                    }
                }
                drawing.pen_drew(drawn);
            } else if old_pos.is_some() {
                needs_save = true;
                drawing.handle_modified_bitmap(&mut frame_images, &mut frame_textures);
                old_pos = None;
                moving_chunk = MovingChunk::from_mask(Pixels::default());
            }
            drawing.show_cursor();
        } else {
            old_pos = None;
            moving_chunk = MovingChunk::from_mask(Pixels::default());

            let (x, y) = mouse_position();
            let m = Vec2::new(x, y);
            const RADIUS: f32 = 30.0;
            let end = Vec2::new(
                RADIUS * 0.65 / 2.0_f32.sqrt(),
                RADIUS * 0.65 * (1.0 + 1.0 / 2.0_f32.sqrt()),
            );
            draw_line(
                m.x + 0.5 * end.x,
                m.y + 0.5 * end.y,
                m.x + 1.05 * end.x,
                m.y + 1.05 * end.y,
                8.0,
                BLACK,
            );
            draw_triangle(
                m,
                m + Vec2::new(0.0, RADIUS),
                m + Vec2::new(RADIUS / 2.0_f32.sqrt(), RADIUS / 2.0_f32.sqrt()),
                WHITE,
            );
            draw_triangle_lines(
                m,
                m + Vec2::new(0.0, RADIUS),
                m + Vec2::new(RADIUS / 2.0_f32.sqrt(), RADIUS / 2.0_f32.sqrt()),
                2.0,
                BLACK,
            );
            draw_line(
                m.x + 0.5 * end.x,
                m.y + 0.5 * end.y,
                m.x + end.x,
                m.y + end.y,
                4.0,
                WHITE,
            );
        }

        next_frame().await
    }
}

#[derive(Clone)]
struct MovingChunk {
    the_mask: Pixels,
    the_chunks: Vec<Pixels>,
    empty: Pixels,
}

impl MovingChunk {
    fn from_mask(mask: Pixels) -> Self {
        MovingChunk {
            the_mask: mask,
            the_chunks: Vec::new(),
            empty: Pixels::default(),
        }
    }
    fn shift_by(&mut self, offset: isize) {
        self.the_mask.shift_by(offset);
        for chunk in self.the_chunks.iter_mut() {
            chunk.shift_by(offset);
        }
    }
    fn insert_chunk(&mut self, which: usize, chunk: Pixels) {
        let chunk: Pixels = chunk.iter().collect();
        while self.the_chunks.len() <= which {
            self.the_chunks.push(Pixels::default());
        }
        self.the_chunks[which] = chunk;
    }
    fn mask(&self) -> &Pixels {
        &self.the_mask
    }
    fn chunks(&self, which: usize) -> &Pixels {
        if let Some(chunk) = self.the_chunks.get(which) {
            &chunk
        } else {
            &self.empty
        }
    }
}
