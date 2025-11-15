//! A module for parsing simple grid data definitions and rendering them as PNG or SVG images.
//! It supports mapping color codes (u8) to specific RGB values and defining rendering parameters.

use image::{Rgb, RgbImage};
use std::collections::HashMap;
use std::error::Error;
use std::fmt::Write;
use std::fs;
use std::path::Path;

/// Represents a pixel grid configuration, containing color mappings, the grid layout,
/// cell size, and grid line visibility settings.
pub struct PixGrid {
    /// Maps the color code (u8) found in `grid_data` to a specific RGB color.
    pub(crate) color_map: HashMap<u8, Rgb<u8>>,
    /// The actual grid data, where each inner Vec<u8> is a row of color codes.
    pub(crate) grid_data: Vec<Vec<u8>>,
    /// The size (width and height) in pixels for each cell in the grid.
    pub(crate) cell_size: u32,
    /// Optional RGB color for drawing grid lines. None means no lines are drawn.
    pub(crate) grid_color: Option<Rgb<u8>>,
}

impl PixGrid {
    /// Parses a string of input data into a PixGrid configuration.
    ///
    /// The input format expects configuration parameters (cell_size, grid_color) and
    /// color definitions (code=R,G,B) followed by a separator "---" and the grid data.
    ///
    /// # Errors
    /// Returns an error if no color definitions are found.
    pub fn parse(input_data: &str) -> Result<PixGrid, Box<dyn Error>> {
        let mut color_map = HashMap::new();
        let mut grid_data = Vec::new();
        let mut reading_grid = false;

        // Default values for parameters
        let mut cell_size: u32 = 10;
        let mut grid_color: Option<Rgb<u8>> = None;

        for line in input_data.lines() {
            // Clean up the line by removing content after '#' (comments) and trimming whitespace.
            let clean_line = if let Some((before_hash, _)) = line.split_once('#') {
                before_hash.trim()
            } else {
                line.trim()
            };

            if clean_line.is_empty() {
                continue;
            }

            // Check for the grid data separator
            if clean_line == "---" {
                reading_grid = true;
                continue;
            }

            if !reading_grid {
                // Handles parameters AND colors (key=value format)
                if let Some((key_str, value_str)) = clean_line.split_once('=') {
                    let key = key_str.trim();
                    let value = value_str.trim();

                    match key {
                        // Handling rendering parameters
                        "cell_size" => {
                            if let Ok(size) = value.parse::<u32>() {
                                cell_size = size;
                            } else {
                                eprintln!("Warning: Invalid 'cell_size': {}", value);
                            }
                        }
                        // New: Handling grid_color parameter (R, G, B)
                        "grid_color" => {
                            let rgb_values: Vec<u8> = value
                                .split(',')
                                .filter_map(|s| s.trim().parse::<u8>().ok())
                                .collect();

                            if rgb_values.len() == 3 {
                                grid_color =
                                    Some(Rgb([rgb_values[0], rgb_values[1], rgb_values[2]]));
                            } else {
                                eprintln!(
                                    "Warning: Invalid 'grid_color' (expected R,G,B): {}",
                                    value
                                );
                            }
                        }
                        "draw_grid" => {
                            eprintln!(
                                "Warning: 'draw_grid' is deprecated. Use 'grid_color = R, G, B' instead."
                            );
                        }

                        // Default case: color handling (codes = R,G,B)
                        _ => {
                            // Try to parse the key as a color code (u8)
                            if let Ok(code) = key.parse::<u8>() {
                                // Parse the value string (e.g., "255,100,0") into RGB components
                                let rgb_values: Vec<u8> = value
                                    .split(',')
                                    .filter_map(|s| s.trim().parse::<u8>().ok())
                                    .collect();

                                if rgb_values.len() == 3 {
                                    let color = Rgb([rgb_values[0], rgb_values[1], rgb_values[2]]);
                                    color_map.insert(code, color);
                                } else {
                                    eprintln!("Warning: Incorrect RGB format for code {}", code);
                                }
                            } else {
                                eprintln!(
                                    "Warning: Unrecognized parameter/color line: {}",
                                    clean_line
                                );
                            }
                        }
                    }
                }
            } else {
                // Reading the grid (space-separated u8 values)
                let row: Vec<u8> = clean_line
                    .split_whitespace()
                    .filter_map(|s| s.parse::<u8>().ok())
                    .collect();

                if !row.is_empty() {
                    grid_data.push(row);
                }
            }
        }

        if color_map.is_empty() {
            return Err("No color definitions found.".into());
        }

        // Returns the complete configuration struct
        Ok(PixGrid {
            color_map,
            grid_data,
            cell_size,
            grid_color,
        })
    }

    /// Generates a PNG image from the stored grid data.
    ///
    /// The image is created by scaling each u8 code in `grid_data` by `cell_size`.
    /// Unknown color codes default to magenta (Rgb([255, 0, 255])).
    ///
    /// # Errors
    /// Returns an error if the grid data is empty or if image saving fails.
    pub fn generate_png(&self, output_path: &Path) -> Result<(), Box<dyn Error>> {
        if self.grid_data.is_empty() || self.grid_data[0].is_empty() {
            return Err("The grid data is empty or invalid.".into());
        }

        let grid_height = self.grid_data.len() as u32;
        let grid_width = self.grid_data[0].len() as u32;
        let cell_size = self.cell_size;

        let image_width = grid_width * cell_size;
        let image_height = grid_height * cell_size;

        let mut img: RgbImage = RgbImage::new(image_width, image_height);
        let unknown_color = Rgb([255, 0, 255]); // Default color for missing codes

        // 1. Iterate through the grid and draw colored cells
        for (y_grid, row) in self.grid_data.iter().enumerate() {
            for (x_grid, &color_code) in row.iter().enumerate() {
                // Get the color, defaulting to unknown_color if the code is not in the map
                let color = self.color_map.get(&color_code).unwrap_or(&unknown_color);

                // Draw a cell of size cell_size x cell_size
                for i in 0..cell_size {
                    for j in 0..cell_size {
                        let x = (x_grid as u32) * cell_size + i;
                        let y = (y_grid as u32) * cell_size + j;

                        img.put_pixel(x, y, *color);
                    }
                }
            }
        }

        // 2. Draw grid lines if the setting is enabled
        if let Some(grid_color) = self.grid_color {
            // Draw vertical lines
            for x_grid_line in 1..grid_width {
                let x_start = x_grid_line * cell_size;
                for y_pixel in 0..image_height {
                    img.put_pixel(x_start, y_pixel, grid_color);
                }
            }

            // Draw horizontal lines
            for y_grid_line in 1..grid_height {
                let y_start = y_grid_line * cell_size;
                for x_pixel in 0..image_width {
                    img.put_pixel(x_pixel, y_start, grid_color);
                }
            }
        }

        println!("Saving the PNG image to: {}", output_path.display());
        img.save(output_path)?;

        Ok(())
    }

    /// Generates an SVG image (XML format) from the stored grid data.
    ///
    /// Each grid cell is drawn as an SVG `<rect>` element.
    ///
    /// # Errors
    /// Returns an error if the grid data is empty or if file writing fails.
    pub fn generate_svg(&self, output_path: &Path) -> Result<(), Box<dyn Error>> {
        if self.grid_data.is_empty() || self.grid_data[0].is_empty() {
            // Corrected error message to English
            return Err("The grid data is empty or invalid.".into());
        }

        let grid_height = self.grid_data.len() as u32;
        let grid_width = self.grid_data[0].len() as u32;
        let cell_size = self.cell_size;

        let image_width = grid_width * cell_size;
        let image_height = grid_height * cell_size;

        let unknown_color = Rgb([255, 0, 255]); // Magenta fallback

        let mut svg_content = String::new();

        // 1. SVG Header: Defines the canvas size and ensures crisp edges for scaling
        writeln!(
            svg_content,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" shape-rendering=\"crispEdges\">",
            image_width, image_height
        )?;

        // 2. Drawing cells (rectangles)
        for (y_grid, row) in self.grid_data.iter().enumerate() {
            for (x_grid, &color_code) in row.iter().enumerate() {
                let color = self.color_map.get(&color_code).unwrap_or(&unknown_color);
                let Rgb([r, g, b]) = *color;
                let fill_color = format!("rgb({},{},{})", r, g, b);

                let x = (x_grid as u32) * cell_size;
                let y = (y_grid as u32) * cell_size;

                // Write the <rect> element for the current cell
                writeln!(
                    svg_content,
                    " 	<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" />",
                    x, y, cell_size, cell_size, fill_color
                )?;
            }
        }

        // 3. Drawing the grid (lines)
        if let Some(grid_color) = self.grid_color {
            let Rgb([r, g, b]) = grid_color;
            let stroke_color = format!("rgb({},{},{})", r, g, b);

            writeln!(
                svg_content,
                " 	<g fill=\"none\" stroke=\"{}\" stroke-width=\"1\">",
                stroke_color
            )?;

            // Vertical lines
            for x_grid_line in 1..grid_width {
                let x = x_grid_line * cell_size;
                writeln!(
                    svg_content,
                    " 	 	<line x1=\"{}\" y1=\"0\" x2=\"{}\" y2=\"{}\" />",
                    x, x, image_height
                )?;
            }

            // Horizontal lines
            for y_grid_line in 1..grid_height {
                let y = y_grid_line * cell_size;
                writeln!(
                    svg_content,
                    " 	 	<line x1=\"0\" y1=\"{}\" x2=\"{}\" y2=\"{}\" />",
                    y, image_width, y
                )?;
            }

            writeln!(svg_content, " 	</g>")?;
        }

        // 4. Closing the SVG
        writeln!(svg_content, "</svg>")?;

        // 5. Writing the file
        println!("Saving the SVG image to: {}", output_path.display());
        fs::write(output_path, svg_content)?;

        Ok(())
    }

    /// Creates a new PixGrid structure by reading a PNG file.
    ///
    /// It samples the top-left pixel of each cell (defined by `cell_size`) to determine
    /// the cell's color and assigns a unique `u8` code to each distinct color found.
    ///
    /// # Arguments
    /// * `input_path` - Path to the PNG file to read.
    /// * `cell_size` - The expected cell size (N x N) in pixels.
    /// * `grid_color_setting` - Sets the `grid_color` parameter for the resulting PixGrid.
    ///
    /// # Errors
    /// Returns an error if the image dimensions are not multiples of `cell_size`,
    /// if the image is empty, or if no colors are found.
    pub fn from_png(
        input_path: &Path,
        cell_size: u32,
        grid_color_setting: Option<Rgb<u8>>,
    ) -> Result<PixGrid, Box<dyn Error>> {
        // 1. Open and prepare the image
        let img = image::open(input_path)?.to_rgb8();
        let width = img.width();
        let height = img.height();

        // Dimension checking
        if width == 0 || height == 0 || cell_size == 0 {
            return Err("Image dimensions or cell_size must be positive.".into());
        }
        // Ensure image dimensions are a clean multiple of cell_size
        if width % cell_size != 0 || height % cell_size != 0 {
            return Err(format!(
                "Image dimensions ({}x{}) must be a multiple of cell_size ({})",
                width, height, cell_size
            )
            .into());
        }

        let grid_width = width / cell_size;
        let grid_height = height / cell_size;

        let mut color_map: HashMap<u8, Rgb<u8>> = HashMap::new();
        // Reverse map to quickly find the code for an existing color
        let mut color_to_code: HashMap<Rgb<u8>, u8> = HashMap::new();
        let mut grid_data: Vec<Vec<u8>> = Vec::with_capacity(grid_height as usize);

        let mut next_code: u8 = 1; // Start codes at 1 to avoid default 0

        // 2. Iterate through the grid, extract colors, and create mappings
        for y_grid in 0..grid_height {
            let mut row: Vec<u8> = Vec::with_capacity(grid_width as usize);

            for x_grid in 0..grid_width {
                // Retrieves the color of the top-left pixel of the cell (sampling)
                let x_pixel = x_grid * cell_size;
                let y_pixel = y_grid * cell_size;

                // FIX: Since `img` is already RgbImage (u8), `get_pixel` returns &Rgb<u8>.
                // We clone the *content* of the reference to get an owned Rgb<u8> value.
                let pixel_color: Rgb<u8> = *img.get_pixel(x_pixel, y_pixel);

                // Note: This implementation assumes all pixels within the cell are the same color.
                // It only checks the top-left pixel.

                // 3. Map the color to a new u8 code (or retrieve existing code)
                let code = *color_to_code.entry(pixel_color).or_insert_with(|| {
                    // It's a new color, assign the next code and update maps
                    let new_code = next_code;
                    next_code = next_code.checked_add(1).unwrap_or_else(|| {
                        eprintln!("Warning: Reached maximum color code (255). Reusing code 255.");
                        255
                    });

                    // Add to the color-to-code and code-to-color maps
                    color_map.insert(new_code, pixel_color);
                    new_code
                });

                row.push(code);
            }
            grid_data.push(row);
        }

        // 4. Building the PixGrid structure
        if color_map.is_empty() {
            return Err("No colors found in PNG to create PixGrid.".into());
        }

        Ok(PixGrid {
            color_map,
            grid_data,
            cell_size,
            grid_color: grid_color_setting,
        })
    }
}

// --- UNIT TESTS ---
#[cfg(test)]
mod tests {
    use super::*;
    use image::Rgb;
    use std::fs;
    use std::path::PathBuf;

    /// Helper function to parse the standard input for tests
    fn get_test_pixgrid() -> PixGrid {
        let input = r#"
# --- RENDERING PARAMETERS ---
cell_size = 40
grid_color = 0, 0, 0 # New: Black grid lines (RGB)

# --- COLOR DEFINITIONS (Code = R, G, B) ---
0 = 255, 255, 255 # White (background)
1 = 0, 0, 0       # Black (border/fill)
2 = 255, 165, 0   # Orange
3 = 0, 255, 0     # Green

# --- GRID START (space-separated) ---
---
0 0 0 0 0 0 0 0 0 0
0 1 1 1 1 1 1 1 1 0
0 1 2 2 2 2 2 2 1 0
0 1 2 3 3 3 3 2 1 0
0 1 2 3 0 0 3 2 1 0
0 1 2 3 0 0 3 2 1 0
0 1 2 3 3 3 3 2 1 0
0 1 2 2 2 2 2 2 1 0
0 1 1 1 1 1 1 1 1 0
0 0 0 0 0 0 0 0 0 0
        "#;
        PixGrid::parse(input).expect("Helper failed to parse standard input")
    }

    /// Helper function to parse input with no grid lines
    fn get_test_pixgrid_no_grid() -> PixGrid {
        let input = r#"
# --- RENDERING PARAMETERS ---
cell_size = 40
# grid_color is omitted

# --- COLOR DEFINITIONS (Code = R, G, B) ---
0 = 255, 255, 255 # White (background)
1 = 0, 0, 0       # Black (border/fill)

# --- GRID START (space-separated) ---
---
0 1 0
1 0 1
0 1 0
        "#;
        PixGrid::parse(input).expect("Helper failed to parse input without grid color")
    }

    #[test]
    /// Tests the PixGrid::parse function with a complete and valid input string.
    fn test_parse_valid_input() {
        let pix_grid = get_test_pixgrid();

        // 1. Test configuration parameters
        assert_eq!(pix_grid.cell_size, 40, "cell_size should be 40");
        // Check grid_color is Some(Black)
        assert_eq!(
            pix_grid.grid_color.unwrap(),
            Rgb([0, 0, 0]),
            "grid_color should be Black (0, 0, 0)"
        );

        // 2. Test grid dimensions
        assert_eq!(pix_grid.grid_data.len(), 10, "Grid should have 10 rows");

        // ... (rest of the assertions)
        assert_eq!(
            pix_grid.color_map.len(),
            4,
            "Color map should have 4 entries"
        );
        assert_eq!(
            *pix_grid.color_map.get(&0).unwrap(),
            Rgb([255, 255, 255]),
            "Color 0 should be White"
        );

        // 4. Test grid data content (sample checks)
        // Check the center cell (4, 4) which should be 0 (White)
        assert_eq!(
            pix_grid.grid_data[4][4], 0,
            "Cell (4, 4) should be color code 0"
        );
    }

    #[test]
    /// Tests parsing when grid_color is omitted (should be None)
    fn test_parse_omitted_grid_color() {
        let pix_grid = get_test_pixgrid_no_grid();
        assert!(
            pix_grid.grid_color.is_none(),
            "grid_color should be None when omitted"
        );
    }

    #[test]
    /// Tests that parsing fails when no color definitions are provided.
    fn test_parse_no_colors() {
        let input = "cell_size = 10\n--- \n0 0 0";
        assert!(
            PixGrid::parse(input).is_err(),
            "Parsing should fail if no colors are defined"
        );
    }

    #[test]
    /// Tests that PNG generation creates a file with the correct dimensions.
    fn test_generate_png_output() {
        let pix_grid = get_test_pixgrid();
        let output_path = PathBuf::from("temp_test_output.png");

        // Ensure cleanup of previous runs and defer cleanup for this run
        let _ = fs::remove_file(&output_path);

        // 1. Generate the PNG file
        pix_grid
            .generate_png(&output_path)
            .expect("PNG generation failed");

        // 2. Assert file was created and has correct size (10x10 grid * 40 cell_size = 400x400)
        assert!(
            output_path.exists(),
            "PNG file should exist after generation"
        );

        match image::open(&output_path) {
            Ok(img) => {
                assert_eq!(img.width(), 400, "PNG width should be 400 pixels");
                assert_eq!(img.height(), 400, "PNG height should be 400 pixels");
            }
            Err(e) => panic!("Failed to open generated PNG for verification: {}", e),
        }

        // 3. Cleanup
        fs::remove_file(&output_path).expect("Failed to clean up generated PNG file");
    }

    #[test]
    /// Tests that SVG generation creates a file with the correct header and stroke color.
    fn test_generate_svg_output() {
        let pix_grid = get_test_pixgrid();
        let output_path = PathBuf::from("temp_test_output.svg");

        // Ensure cleanup of previous runs and defer cleanup for this run
        let _ = fs::remove_file(&output_path);

        // 1. Generate the SVG file
        pix_grid
            .generate_svg(&output_path)
            .expect("SVG generation failed");

        // 2. Assert file was created
        assert!(
            output_path.exists(),
            "SVG file should exist after generation"
        );

        // 3. Assert content (400x400) and the stroke color
        let content = fs::read_to_string(&output_path).expect("Failed to read generated SVG file");

        let expected_header = r#"<svg xmlns="http://www.w3.org/2000/svg" width="400" height="400" shape-rendering="crispEdges">"#;
        let expected_stroke = r#"stroke="rgb(0,0,0)""#;

        assert!(
            content.contains(expected_header),
            "SVG header should contain correct dimensions (400x400)"
        );
        assert!(
            content.contains(expected_stroke),
            "SVG stroke color should be rgb(0,0,0) based on input grid_color"
        );

        // 4. Cleanup
        fs::remove_file(&output_path).expect("Failed to clean up generated SVG file");
    }
}
