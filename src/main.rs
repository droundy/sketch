use std::{f32::consts::PI, sync::atomic::AtomicBool};

use macroquad::prelude::{
    draw_circle, draw_circle_lines, draw_line, draw_rectangle, draw_rectangle_lines, draw_texture,
    is_key_pressed, is_mouse_button_down, is_mouse_button_pressed, mouse_position, next_frame,
    screen_height, screen_width, Color, Conf, KeyCode, MouseButton, Vec2, BLACK, GRAY, WHITE,
};
use macroquad::shapes::draw_poly;
use macroquad::ui::root_ui;

mod layer;
mod tween;
use layer::Layer;

fn conf() -> Conf {
    Conf {
        window_title: String::from("Fun draw"),
        window_width: 1920,
        window_height: 1080,
        fullscreen: true,
        ..Default::default()
    }
}

fn color_selector_color(fx: f32, fy: f32) -> Option<Color> {
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
        let rgb = rgb.map(|c| c * (1.0 - extra) + extra);
        Some(Color {
            r: rgb[0],
            g: rgb[1],
            b: rgb[2],
            a: 1.0,
        })
    } else if radius <= 1.0 {
        let x = 1.0 - (radius - RSATURATED) / (1.0 - RSATURATED);
        let r = x * (1.0 - (1.0 - x).powi(2));
        Some(Color {
            r: r * rgb[0],
            g: r * rgb[1],
            b: r * rgb[2],
            a: 1.0,
        })
    } else {
        None
    }
}

fn color_selector(color: &mut Color) -> bool {
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
                draw_rectangle(x, y, dx, dy, c);
            }
        }
    }

    static AM_DRAGGING: AtomicBool = AtomicBool::new(false);
    if is_mouse_button_down(MouseButton::Left) {
        let pos = mouse_position();
        let fx = (pos.0 + w - swidth) / w;
        let fy = (pos.1 + h - sheight) / h;
        if let Some(c) = color_selector_color(fx, fy) {
            *color = c;
            AM_DRAGGING.store(true, std::sync::atomic::Ordering::Relaxed);
            return true;
        } else if AM_DRAGGING.load(std::sync::atomic::Ordering::Relaxed) {
            return true;
        }
    } else {
        AM_DRAGGING.store(false, std::sync::atomic::Ordering::Relaxed);
    }
    false
}

impl Drawing {
    fn frame_selector(&mut self) -> bool {
        self.layers[self.current].frame_selector(&mut self.time)
    }
    fn get_frame_data_mut(&mut self) -> &mut [[u8; 4]] {
        self.layers[self.current].get_frame_data_mut(self.time)
    }
    fn handle_modified_bitmap(&mut self) {
        self.layers[self.current].handle_modified_bitmap(self.time);
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
            if !self.am_animating {
                self.time = self.layers[self.current].closest_time(self.time);
            }
            true
        } else {
            false
        }
    }
    fn layer_selector(&mut self) -> bool {
        const WIDTH: f32 = 50.0;
        const HEIGHT: f32 = 40.0;
        const NUM_TOOLS: usize = 3;
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
            let y = (NUM_TOOLS + 1 + i) as f32 * HEIGHT;
            draw_rectangle(0.0, y, WIDTH, HEIGHT, l.color);
            outline(y, i == self.current);
        }
        let y = (NUM_TOOLS + 1 + self.layers.len()) as f32 * HEIGHT;
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
        let y = 0.0;
        outline(y, self.tool == Tool::BigPen);
        draw_circle(WIDTH * 0.5, y + HEIGHT * 0.5, 10.0, WHITE);
        outline(y + HEIGHT, self.tool == Tool::LittlePen);
        draw_circle(WIDTH * 0.5, y + HEIGHT * 1.5, 2.0, WHITE);
        outline(y + 2.0 * HEIGHT, self.tool == Tool::Eraser);
        draw_circle_lines(WIDTH * 0.5, y + HEIGHT * 2.5, 10.0, 2.0, WHITE);
        if is_mouse_button_down(MouseButton::Left) || is_mouse_button_pressed(MouseButton::Left) {
            let (x, y) = mouse_position();
            let y = (y / HEIGHT) as usize;
            if x < WIDTH {
                if is_mouse_button_pressed(MouseButton::Left) {
                    if y == NUM_TOOLS + 1 + self.layers.len() {
                        self.layers.push(Layer::new(self.time));
                        self.current = self.layers.len() - 1;
                        if self.tool == Tool::Eraser {
                            self.tool = Tool::BigPen;
                        }
                    } else if y > NUM_TOOLS && y < NUM_TOOLS + 1 + self.layers.len() {
                        self.current = y - NUM_TOOLS - 1;
                        if self.tool == Tool::Eraser {
                            self.tool = Tool::BigPen;
                        }
                    } else if y == 0 {
                        self.tool = Tool::BigPen;
                    } else if y == 1 {
                        self.tool = Tool::LittlePen;
                    } else if y == 2 {
                        self.tool = Tool::Eraser;
                    }
                }
                return true;
            }
        }
        false
    }
}

#[derive(PartialEq, Eq, Debug, Copy, Clone)]
enum Tool {
    Eraser,
    BigPen,
    LittlePen,
}

struct Drawing {
    am_animating: bool,
    current: usize,
    height: u16,
    width: u16,
    time: f32,
    tool: Tool,
    layers: Vec<Layer>,
}

#[macroquad::main(conf)]
async fn main() {
    let mut old_pos: Option<Vec2> = None;
    let mut drawing = Drawing {
        am_animating: false,
        time: 0.0,
        current: 0,
        tool: Tool::BigPen,
        width: screen_width() as u16,
        height: screen_height() as u16,
        layers: vec![Layer::new(0.0)],
    };
    let width = drawing.width as usize;
    let height = drawing.height as usize;
    loop {
        // clear_background(WHITE);
        if is_key_pressed(KeyCode::Escape) {
            return;
        }
        if drawing.am_animating {
            drawing.time += 0.01;
            if drawing.time > 1.0 {
                drawing.time = 0.0;
            }
        }

        for l in drawing.layers.iter_mut() {
            draw_texture(l.texture(drawing.time), 0.0, 0.0, l.color);
        }
        let animation_button_selected = drawing.animation_button();
        let color_selected = color_selector(&mut drawing.layers[drawing.current].color);
        let frame_selected = drawing.frame_selector();
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
                let radius = match drawing.tool {
                    Tool::LittlePen => 2.0,
                    Tool::BigPen => 10.0,
                    Tool::Eraser => 20.0,
                };
                let drawn = if drawing.tool == Tool::Eraser {
                    [0; 4]
                } else {
                    [255; 4]
                };
                let data = drawing.get_frame_data_mut();
                if let Some(old) = old_pos {
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
                                data[x + y * width] = drawn;
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
                            data[x + y * width] = drawn;
                        }
                        old_pos = Some(pos);
                    }
                }
                drawing.handle_modified_bitmap();
            } else {
                old_pos = None;
            }
        } else {
            old_pos = None;
        }

        next_frame().await
    }
}
