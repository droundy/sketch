use macroquad::{prelude::*, ui::root_ui};

fn conf() -> Conf {
    Conf {
        window_title: String::from("Fun draw"),
        window_width: 1260,
        window_height: 768,
        fullscreen: true,
        ..Default::default()
    }
}

fn color_selector(color: &mut Color) {
    if root_ui().button(None, "White") {
        *color = WHITE;
    }
    if root_ui().button(None, "Red") {
        *color = RED;
    }
    if root_ui().button(None, "Blue") {
        *color = BLUE;
    }
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

        color_selector(&mut color);
        root_ui().slider(0, "time", 0.0..1.0, &mut time);

        if !root_ui().is_mouse_captured() {
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
        if is_key_pressed(KeyCode::Escape) {
            return;
        }

        draw_texture(Texture2D::from_image(&image), 0.0, 0.0, color);
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
