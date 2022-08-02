use std::f32::consts::PI;

use macroquad::{
    prelude::*,
    ui::{root_ui, InputHandler},
};

fn conf() -> Conf {
    Conf {
        window_title: String::from("Fun draw"),
        window_width: 1260,
        window_height: 768,
        fullscreen: true,
        ..Default::default()
    }
}

fn color_selector_color(fx: f32, fy: f32) -> Color {
    let fx = fx * 2.0 - 1.0;
    let fy = fy * 2.0 - 1.0;
    let angle = fy.atan2(fx) + PI;
    let radius = (fx * fx + fy * fy).sqrt();
    const PI6: f32 = PI / 3.0;
    let c1 = (((angle + PI6) % (2.0 * PI / 3.0)) / (PI / 3.0) - 1.0).abs();
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
    let radius = if radius > 0.8 {
        (2.0_f32).sqrt() * radius.sqrt() / (1.0 + c1 * c1).sqrt()
    } else {
        (2.0_f32).sqrt() * radius.sqrt() / (1.0 + c1 * c1).sqrt()
    
    };
    Color {
        r: radius * rgb[0],
        g: radius * rgb[1],
        b: radius * rgb[2],
        a: 1.0,
    }
}

fn color_selector(color: &mut Color) -> bool {
    let swidth = screen_width();
    let sheight = screen_height();

    let radius = 100.0;
    let cx = swidth - radius;
    let cy = sheight - radius;
    for theta in (0..100).map(|i| i as f32 * (PI / 50.0)) {
        let (sin, cos) = theta.sin_cos();
        draw_line(
            cx,
            cy,
            cx + radius * cos,
            cy + radius * sin,
            radius / 100.0,
            color_selector_color(sin, cos),
        );
    }

    let w = 200.0;
    let h = 200.0;
    draw_rectangle(swidth - w, sheight - h, w, h, RED);
    let dx = w * 0.01;
    let dy = h * 0.01;
    for i in 0..100 {
        for j in 0..100 {
            let fx = i as f32 * 0.01;
            let fy = j as f32 * 0.01;
            let x = fx * w + swidth - w;
            let y = fy * h + sheight - h;
            draw_rectangle(x, y, dx, dy, color_selector_color(fx, fy));
        }
    }

    // let radius = 100.0;
    // let center = Vec2::new(swidth - radius, sheight - radius);
    // const DTHETA: f32 = PI / 50.0;
    // for theta in (0..100).map(|i| i as f32 * DTHETA) {
    //     for r in (0..100).rev().map(|i| i as f32 * radius * 0.01) {
    //         let (sin, cos) = theta.sin_cos();
    //         let v1 = center + r * Vec2::new(sin, cos);
    //         let (nsin, ncos) = (theta + DTHETA).sin_cos();
    //         let v2 = center + r * Vec2::new(nsin, ncos);
    //         draw_triangle(v1, v2, center, color_selector_color(sin, cos));
    //     }
    // }
    if is_mouse_button_down(MouseButton::Left) {
        let pos = mouse_position();
        if pos.0 > swidth - w && pos.1 > sheight - h {
            let fx = (pos.0 + w - swidth) / w;
            let fy = (pos.1 + h - sheight) / h;
            *color = color_selector_color(fx, fy);
            return true;
        }
    }
    false
}
#[macroquad::main(conf)]
async fn main() {
    let mut points = Vec::new();
    let mut old_pos: Option<Vec2> = None;
    let mut time = 0.0;
    let mut image = Image::gen_image_color(2000, 1024, BLACK);
    let width = image.width as usize;
    let radius = 5.0;
    let mut color = WHITE;
    loop {
        // clear_background(WHITE);
        if is_key_pressed(KeyCode::Escape) {
            return;
        }

        root_ui().slider(0, "time", 0.0..1.0, &mut time);

        draw_texture(Texture2D::from_image(&image), 0.0, 0.0, color);
        if !root_ui().is_mouse_captured() && !color_selector(&mut color) {
            if is_mouse_button_down(MouseButton::Left) {
                let pos = mouse_position();
                let pos = Vec2::new(pos.0, pos.1);
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
                        image.height as usize - 1 - (radius as usize),
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
                                image.get_image_data_mut()[x + y * width] = [255, 255, 255, 255];
                            }
                        }
                    }
                }
                let x_start = std::cmp::max(radius as usize, pos.x as usize) - radius as usize;
                let x_stop = std::cmp::min(width - 1 - (radius as usize), pos.x as usize)
                    + radius as usize
                    + 1;
                let y_start = std::cmp::max(radius as usize, pos.y as usize) - radius as usize;
                let y_stop = std::cmp::min(
                    image.height as usize - 1 - (radius as usize),
                    pos.y as usize,
                ) + radius as usize
                    + 1;
                for x in x_start..x_stop {
                    for y in y_start..y_stop {
                        if (x as f32 - pos.x).powi(2) + (y as f32 - pos.y).powi(2) < radius.powi(2)
                        {
                            image.get_image_data_mut()[x + y * width] = [255, 255, 255, 255];
                        }
                        old_pos = Some(pos);
                    }
                }
                points.push(pos);
            } else {
                old_pos = None;
            }
        }

        // for p in points.iter() {
        //     draw_circle(p.0, p.1, 5.0, BLUE);
        // }
        // draw_line(40.0, 40.0, 100.0, 200.0, 15.0, BLUE);
        // draw_rectangle(screen_width() / 2.0 - 60.0, 100.0, 120.0, 60.0, GREEN);
        // draw_circle(screen_width() - 30.0, screen_height() - 30.0, 15.0, YELLOW);

        // draw_text("IT WORKS!", 20.0, 20.0, 30.0, DARKGRAY);

        next_frame().await
    }
}
