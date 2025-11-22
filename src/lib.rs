//! A module for parsing simple grid data definitions and rendering them as PNG or SVG images.
//! It supports mapping color codes (u8) to specific RGB values and defining rendering parameters.
//! ```rust
//! use pixgrid::PixGrid;
//! use std::error::Error;
//! use std::fs;
//! use std::path::Path;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     let input_file_path = "instances/all_diff.pg";
//!
//!     // Modify here to choose the output extension and format.
//!     let output_file_name = "nested_squares.svg";
//!
//!     let output_path = Path::new(output_file_name);
//!
//!     // 1. Reading the file
//!     println!("Reading file: {}", input_file_path);
//!     let input_data = fs::read_to_string(input_file_path).map_err(|e| {
//!         format!(
//!             "Error reading file {}: {}",
//!             input_file_path, e
//!         )
//!     })?;
//!
//!     // 2. Parsing colors and the grid (PixGrid contains everything)
//!     // Using the static method PixGrid::parse
//!     let pg = PixGrid::parse(&input_data)?;
//!
//!     println!(
//!         "Defined colors: {:?}",
//!         pg.color_map.keys().collect::<Vec<_>>()
//!     );
//!     println!(
//!         "Grid read: {} rows x {} columns",
//!         pg.grid_data.len(),
//!         pg.grid_data.first().map_or(0, |r| r.len())
//!     );
//!
//!     // 3. Image generation (selection logic)
//!     let extension = output_path
//!         .extension()
//!         .and_then(std::ffi::OsStr::to_str)
//!         .unwrap_or_default();
//!
//!     match extension {
//!         "png" => {
//!             println!("Generating in PNG format...");
//!             // Using the instance method
//!             pg.export_png(output_path)?;
//!         }
//!         "svg" => {
//!             println!("Generating in SVG format...");
//!             // Using the instance method
//!             pg.export_svg(output_path)?;
//!         }
//!         _ => {
//!             return Err(format!(
//!                 "Unsupported file extension: '{}'. Use '.png' or '.svg'.",
//!                 extension
//!             )
//!             .into());
//!         }
//!     }
//!
//!     println!("Generation finished successfully!");
//!     Ok(())
//! }
//! ```

use image::{Rgb, RgbImage};
use rand::rng;
use rand::seq::SliceRandom;
use rand::Rng;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::iter;
use std::path::Path;

/// Utility function to convert HSV to RGB.
/// h: Hue [0.0, 360.0], s: Saturation [0.0, 1.0], v: Value [0.0, 1.0]
fn hsv_to_rgb(h: f32, s: f32, v: f32) -> Rgb<u8> {
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;

    let (r_prime, g_prime, b_prime) = if h < 60.0 {
        (c, x, 0.0)
    } else if h < 120.0 {
        (x, c, 0.0)
    } else if h < 180.0 {
        (0.0, c, x)
    } else if h < 240.0 {
        (0.0, x, c)
    } else if h < 300.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };

    Rgb([
        ((r_prime + m) * 255.0).round() as u8,
        ((g_prime + m) * 255.0).round() as u8,
        ((b_prime + m) * 255.0).round() as u8,
    ])
}
fn color_distance(c1: &Rgb<u8>, c2: &Rgb<u8>) -> f32 {
    let r_diff = c1[0] as f32 - c2[0] as f32;
    let g_diff = c1[1] as f32 - c2[1] as f32;
    let b_diff = c1[2] as f32 - c2[2] as f32;
    (r_diff * r_diff + g_diff * g_diff + b_diff * b_diff).sqrt()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CellPosition {
    pub x: usize, // Index de Colonne (Largeur)
    pub y: usize, // Index de Ligne (Hauteur)
}

/// Represents a pixel grid configuration, containing color mappings, the grid layout,
/// cell size, and grid line visibility settings.
#[derive(Debug, Clone)]
pub struct PixGrid {
    /// Maps the color code (u8) found in `grid_data` to a specific RGB color.
    pub color_map: HashMap<u8, Rgb<u8>>,
    /// The actual grid data, where each inner Vec<u8> is a row of color codes.
    pub grid_data: Vec<Vec<u8>>,
    /// The size (width and height) in pixels for each cell in the grid.
    pub cell_size: u32,
    /// Optional RGB color for drawing grid lines. None means no lines are drawn.
    pub grid_color: Option<Rgb<u8>>,
}

impl PixGrid {
    /// ✨ **Constructor:** Creates a new `PixGrid` with the specified dimensions.
    ///
    /// # Arguments
    ///
    /// * `width` - The number of columns (grid width).
    /// * `height` - The number of rows (grid height).
    ///
    /// # Default Values Used
    ///
    /// * **`cell_size`**: 1 (Default pixel size for each cell).
    /// * **`grid_color`**: `None` (No grid lines drawn by default).
    /// * **`default_color_code`**: 0 (Default color code used to fill the grid).
    pub fn new(width: usize, height: usize) -> Self {
        // Define the default color code for initial filling (e.g., 0 for background/empty)
        const DEFAULT_COLOR_CODE: u8 = 0;

        // 1. Create a row (Vec<u8>) of 'width' size, filled with the default color code.
        let row = vec![DEFAULT_COLOR_CODE; width];

        // 2. Create the vector of vectors (grid_data) of 'height',
        //    where each element is a copy of the 'row'.
        let grid_data = vec![row; height];

        PixGrid {
            color_map: HashMap::new(), // Empty color palette by default
            grid_data,
            cell_size: 1,     // Default cell size of 1 pixel
            grid_color: None, // No grid lines by default
        }
    }

    /// Sets the color code of a single cell at the given logical position.
    ///
    /// Returns Ok(()) on success or Err if the position is out of bounds.
    pub fn set_cell(&mut self, position: CellPosition, color_code: u8) -> Result<(), &'static str> {
        let (width, height) = self.dimensions();

        // Check if the coordinates are within the grid bounds
        if position.y >= height || position.x >= width {
            return Err("Cell position is out of bounds.");
        }

        // Access the grid using the named fields (y then x)
        self.grid_data[position.y][position.x] = color_code;
        Ok(())
    }

    /// Retrieves the color code of a single cell at the given logical position.
    ///
    /// Returns Option<u8>, which is None if the position is out of bounds.
    pub fn get_cell(&self, position: CellPosition) -> Option<u8> {
        // Use the Option trait functions for safe, concise access
        self.grid_data
            .get(position.y) // Checks if the row index (y) is valid
            .and_then(|row| row.get(position.x)) // Checks if the column index (x) is valid
            .copied() // Converts &u8 to u8, returning None if any step failed
    }

    /// Retrieves the position and color code of the 8 surrounding cells (Moore neighborhood)
    /// for the given position. Out-of-bounds cells are skipped.
    ///
    /// Returns a Vec<(CellPosition, u8)>.
    pub fn get_moore_neighbors(&self, center_pos: CellPosition) -> Vec<(CellPosition, u8)> {
        let (width, height) = self.dimensions();
        let mut neighbors = Vec::with_capacity(8);

        // If the center position itself is out of bounds, return an empty list.
        if center_pos.x >= width || center_pos.y >= height {
            return neighbors;
        }

        // Iterate over the 8 directions (dx, dy) around the center
        for dy in -1..=1 {
            for dx in -1..=1 {
                // Skip the central cell (dx=0, dy=0)
                if dx == 0 && dy == 0 {
                    continue;
                }

                // Calculate the potential neighbor's position
                // Use i64 to safely calculate coordinates when subtracting 1 (-1)
                let new_x = center_pos.x as i64 + dx;
                let new_y = center_pos.y as i64 + dy;

                // 1. Check if the coordinates are within the limits [0, width/height - 1]
                if new_x >= 0 && new_x < width as i64 && new_y >= 0 && new_y < height as i64 {
                    let neighbor_pos = CellPosition {
                        x: new_x as usize,
                        y: new_y as usize,
                    };

                    // 2. Retrieve the color code. We can use get_cell for safety.
                    if let Some(code) = self.get_cell(neighbor_pos) {
                        neighbors.push((neighbor_pos, code));
                    }
                }
            }
        }
        neighbors
    }

    /// Retrieves the position and color code of the 4 cardinal neighboring cells
    /// (Von Neumann neighborhood: Top, Bottom, Left, Right) for the given position.
    /// Out-of-bounds cells are skipped.
    ///
    /// Returns a Vec<(CellPosition, u8)>.
    pub fn get_von_neumann_neighbors(&self, center_pos: CellPosition) -> Vec<(CellPosition, u8)> {
        let (width, height) = self.dimensions();
        let mut neighbors = Vec::with_capacity(4);

        // If the center position itself is out of bounds, return an empty list.
        if center_pos.x >= width || center_pos.y >= height {
            return neighbors;
        }

        // Define the 4 cardinal directions: (dx, dy)
        let directions = [
            (0, -1), // Top
            (0, 1),  // Bottom
            (-1, 0), // Left
            (1, 0),  // Right
        ];

        for (dx, dy) in directions.iter() {
            // Use i64 for safe coordinate arithmetic
            let new_x = center_pos.x as i64 + dx;
            let new_y = center_pos.y as i64 + dy;

            // 1. Check if the potential neighbor is within the grid limits
            if new_x >= 0 && new_x < width as i64 && new_y >= 0 && new_y < height as i64 {
                let neighbor_pos = CellPosition {
                    x: new_x as usize,
                    y: new_y as usize,
                };

                // 2. Retrieve the color code and push the neighbor if successful
                if let Some(code) = self.get_cell(neighbor_pos) {
                    neighbors.push((neighbor_pos, code));
                }
            }
        }
        neighbors
    }

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

        // Returns the complete configuration struct
        Ok(PixGrid {
            color_map,
            grid_data,
            cell_size,
            grid_color,
        })
    }

    /// Retrieves the width and height of the grid (in number of cells).
    /// Returns (width, height). Returns (0, 0) if the grid is empty.
    pub fn dimensions(&self) -> (usize, usize) {
        let height = self.grid_data.len();

        let width = if height > 0 {
            // Use the length of the first row for the width.
            // map_or(0, |r| r.len()) safely handles the case of an empty inner vector.
            self.grid_data.first().map_or(0, |r| r.len())
        } else {
            0
        };

        (width, height)
    }

    /// Initializes the grid with a specific color code based on a given density (probability).
    /// Cells not assigned the activation code will remain at the default code (0).
    ///
    /// # Arguments
    /// * `cell_state`: The color code (e.g., 1) to use for the randomized cells.
    /// * `density`: The probability (0.0 to 1.0) that a cell will be set to `cell_state`.
    pub fn random_sparse_fill(&mut self, cell_state: u8, density: f64) {
        let (width, height) = self.dimensions();
        let mut rng = rng();

        // Ensure density is within a valid range
        let p = density.clamp(0.0, 1.0);

        // Iterate through all cells
        for y in 0..height {
            for x in 0..width {
                // Check if a random floating-point number is less than the density
                if rng.random::<f64>() < p {
                    // Set the cell to the activation code
                    self.grid_data[y][x] = cell_state;
                } else {
                    // Set the cell to the default code (usually 0, but explicit reset is safer)
                    self.grid_data[y][x] = 0;
                }
            }
        }
    }

    /// Fills the grid randomly according to a given vector of code proportions.
    ///
    /// # Arguments
    /// * `proportions`: A slice of tuples `(color_code: u8, proportion: f64)`.
    ///   The sum of all proportions must be close to 1.0.
    ///
    /// # Errors
    /// Returns an `Err` if the proportions are invalid (sum is not close to 1.0).
    pub fn random_fill_with_proportions(
        &mut self,
        proportions: &[(u8, f64)],
    ) -> Result<(), &'static str> {
        let (width, height) = self.dimensions();
        let total_cells = width * height;

        if total_cells == 0 {
            return Ok(());
        }

        // --- 1. Validate Proportion Sum ---
        let sum_of_proportions: f64 = proportions.iter().map(|(_, p)| p).sum();
        if (sum_of_proportions - 1.0).abs() > 1e-6 {
            return Err("The sum of proportions must be close to 1.0 (e.g., 100%).");
        }

        // --- 2. Build the List of Codes ---
        let mut codes_list: Vec<u8> = Vec::with_capacity(total_cells);

        for &(code, proportion) in proportions.iter() {
            // Calculate the number of cells required for this code
            let count = (total_cells as f64 * proportion).round() as usize;

            // Add the color code the required number of times using the concise method
            codes_list.extend(iter::repeat_n(code, count));
        }

        // Adjust size if rounding caused a slight difference
        if codes_list.len() > total_cells {
            codes_list.truncate(total_cells);
        } else if codes_list.len() < total_cells {
            // Fill the small difference with the first code (or the most frequent)
            let diff = total_cells - codes_list.len();
            if let Some(&(first_code, _)) = proportions.first() {
                codes_list.extend(std::iter::repeat_n(first_code, diff));
            }
        }

        // --- 3. Random Shuffle ---
        let mut rng = rng();
        codes_list.shuffle(&mut rng);

        // --- 4. Fill the Grid ---
        let mut code_iter = codes_list.into_iter();

        // The inner loop iterates x (columns), the outer loop iterates y (rows)
        for row in self.grid_data.iter_mut() {
            for cell in row.iter_mut() {
                // `next()` should always succeed here as we adjusted the size
                if let Some(code) = code_iter.next() {
                    *cell = code;
                }
            }
        }

        Ok(())
    }

    /// Exports the current PixGrid configuration to a file in the *.pg format.
    ///
    /// The output format includes:
    /// 1. Rendering parameters (cell_size, grid_color).
    /// 2. Color definitions (Code = R, G, B).
    /// 3. The grid data, separated by "---".
    ///
    /// # Arguments
    /// * `output_path` - The path to the file where the data will be written.
    ///
    /// # Returns
    /// A `Result` indicating success or an error if writing to the file fails.
    pub fn export_pg(&self, output_path: &Path) -> Result<(), Box<dyn Error>> {
        // Create or truncate the output file.
        let mut file = File::create(output_path)?;

        // --- 1. Write Rendering Parameters ---
        writeln!(file, "cell_size = {}", self.cell_size)?;

        // Write grid_color if it is defined.
        if let Some(Rgb([r, g, b])) = self.grid_color {
            writeln!(file, "grid_color = {}, {}, {}", r, g, b)?;
        } else {
            // Optional: Write a comment indicating no grid color is set, or simply omit the line.
            writeln!(
                file,
                "# grid_color = <R, G, B> # Omitted: No grid lines are drawn."
            )?;
        }

        if !self.color_map.is_empty() {
            // --- 2. Write Color Definitions ---
            // Sort the color map by code for a consistent and readable output.
            let mut sorted_colors: Vec<(&u8, &Rgb<u8>)> = self.color_map.iter().collect();
            sorted_colors.sort_by_key(|(code, _)| *code);

            for (code, Rgb([r, g, b])) in sorted_colors {
                // We use format! to construct the R, G, B string
                writeln!(file, "{} = {}, {}, {}", code, r, g, b)?;
            }
        }

        // --- 3. Write Grid Data ---
        writeln!(file, "---")?;

        // Write each row of the grid data.
        for row in &self.grid_data {
            // Convert the Vec<u8> row into a space-separated string.
            let row_string: String = row
                .iter()
                .map(|&code| code.to_string())
                .collect::<Vec<String>>()
                .join(" ");

            writeln!(file, "{}", row_string)?;
        }

        // Ensure all buffered content is written to the file.
        file.flush()?;

        Ok(())
    }

    pub fn get_missing_codes(&self) -> Vec<u8> {
        let mut unique_codes = HashSet::new();
        for row in &self.grid_data {
            for &code in row {
                unique_codes.insert(code);
            }
        }
        unique_codes
            .into_iter()
            .filter(|code| !self.color_map.contains_key(code))
            .collect()
    }

    /// Returns a list of standard distinct colors
    fn get_standard_palette() -> Vec<Rgb<u8>> {
        vec![
            Rgb([255, 255, 255]), // White
            Rgb([0, 0, 0]),       // Black
            Rgb([31, 119, 180]),  // Blue
            Rgb([255, 127, 14]),  // Orange
            Rgb([44, 160, 44]),   // Green
            Rgb([214, 39, 40]),   // Red
            Rgb([148, 103, 189]), // Purple
            Rgb([140, 86, 75]),   // Brown
            Rgb([227, 119, 194]), // Pink
            Rgb([127, 127, 127]), // Gray
            Rgb([188, 189, 34]),  // Olive
            Rgb([23, 190, 207]),  // Cyan
            Rgb([255, 215, 0]),   // Gold
            Rgb([0, 255, 127]),   // Spring Green
        ]
    }

    /// Generates a complete color map.
    /// 1. Tries to use standard distinct colors first.
    /// 2. Falls back to Golden Angle approximation for remaining codes.
    fn prepare_complete_color_map(&self) -> HashMap<u8, Rgb<u8>> {
        let mut working_map = self.color_map.clone();
        let mut missing_codes = self.get_missing_codes();

        if missing_codes.is_empty() {
            return working_map;
        }

        // SORTING IS CRITICAL for deterministic results
        missing_codes.sort();

        let palette = Self::get_standard_palette();
        let min_dist = 30.0; // Distance threshold to consider a color "too close"

        // List for codes that couldn't be satisfied by the palette
        let mut codes_for_fallback = Vec::new();

        // --- PHASE 1: Standard Palette ---
        for &code in &missing_codes {
            let mut assigned = false;

            // Try to find a palette color that isn't close to ANY existing color in the map
            for candidate in &palette {
                let mut collision = false;
                for existing in working_map.values() {
                    if color_distance(candidate, existing) < min_dist {
                        collision = true;
                        break;
                    }
                }

                if !collision {
                    working_map.insert(code, *candidate);
                    assigned = true;
                    break;
                }
            }

            if !assigned {
                codes_for_fallback.push(code);
            }
        }

        // --- PHASE 2: Golden Angle Fallback ---
        // If there are still codes without colors, we use your Golden Angle logic.

        if !codes_for_fallback.is_empty() {
            // We iterate strictly on the remaining codes

            // Let's start at 0 (Red) or an arbitrary offset.
            // (Note: Since we might have added Red from the palette, overlap is possible strictly speaking,
            // but Golden Angle quickly diverges).
            let mut current_hue = 0.0;

            // The Golden Angle approx (137.5 degrees).
            let golden_angle = 137.508;

            for &code in &codes_for_fallback {
                // Generate color
                // Saturation 0.8, Value 0.95 = Bright and distinct.
                let generated_color = hsv_to_rgb(current_hue, 0.8, 0.95);

                // We insert directly (assuming the user accepts the mathematical guarantee of distribution)
                // or we could add a check here too, but let's stick to the requested logic logic.
                working_map.insert(code, generated_color);

                // Jump by the golden angle for the next code
                current_hue = (current_hue + golden_angle) % 360.0;
            }
        }

        working_map
    }

    /// Export a PNG image from the stored grid data.
    ///
    /// The image is created by scaling each u8 code in `grid_data` by `cell_size`.
    /// Unknown color codes default to magenta (Rgb([255, 0, 255])).
    ///
    /// # Errors
    /// Returns an error if the grid data is empty or if image saving fails.
    pub fn export_png(&self, output_path: &Path) -> Result<(), Box<dyn Error>> {
        if self.grid_data.is_empty() || self.grid_data[0].is_empty() {
            return Err("The grid data is empty or invalid.".into());
        }

        let complete_map = self.prepare_complete_color_map();

        let grid_height = self.grid_data.len() as u32;
        let grid_width = self.grid_data[0].len() as u32;
        let cell_size = self.cell_size;

        let image_width = grid_width * cell_size;
        let image_height = grid_height * cell_size;

        let mut img: RgbImage = RgbImage::new(image_width, image_height);

        // 1. Iterate through the grid and draw colored cells
        for (y_grid, row) in self.grid_data.iter().enumerate() {
            for (x_grid, &color_code) in row.iter().enumerate() {
                // Get the color, defaulting to unknown_color if the code is not in the map
                let color = complete_map.get(&color_code).unwrap();
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

    /// Export an SVG image (XML format) from the stored grid data.
    ///
    /// Each grid cell is drawn as an SVG `<rect>` element.
    ///
    /// # Errors
    /// Returns an error if the grid data is empty or if file writing fails.
    pub fn export_svg(&self, output_path: &Path) -> Result<(), Box<dyn Error>> {
        if self.grid_data.is_empty() || self.grid_data[0].is_empty() {
            // Corrected error message to English
            return Err("The grid data is empty or invalid.".into());
        }

        let complete_map = self.prepare_complete_color_map();

        let grid_height = self.grid_data.len() as u32;
        let grid_width = self.grid_data[0].len() as u32;
        let cell_size = self.cell_size;

        let image_width = grid_width * cell_size;
        let image_height = grid_height * cell_size;

        let mut svg_content: Vec<u8> = Vec::new(); // Initialisation

        // 1. SVG Header: Defines the canvas size and ensures crisp edges for scaling
        writeln!(
            svg_content,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" shape-rendering=\"crispEdges\">",
            image_width, image_height
        )?;

        // 2. Drawing cells (rectangles)
        for (y_grid, row) in self.grid_data.iter().enumerate() {
            for (x_grid, &color_code) in row.iter().enumerate() {
                let color = complete_map.get(&color_code).unwrap();
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
    use std::io::Read;
    use std::path::PathBuf;
    use std::{env, fs};

    #[test]
    /// Tests the PixGrid::new(width, height) constructor to ensure default values and size are correct.
    fn test_new() {
        const WIDTH: usize = 5;
        const HEIGHT: usize = 3;

        let pg = PixGrid::new(WIDTH, HEIGHT);

        // 1. Test configuration parameters defaults
        assert_eq!(pg.cell_size, 1, "cell_size should default to 1");
        assert!(pg.grid_color.is_none(), "grid_color should default to None");
        assert!(
            pg.color_map.is_empty(),
            "color_map should be empty by default"
        );

        // 2. Test grid dimensions
        assert_eq!(
            pg.grid_data.len(),
            HEIGHT,
            "Grid should have the specified height ({} rows)",
            HEIGHT
        );
        assert!(
            pg.grid_data.iter().all(|row| row.len() == WIDTH),
            "All rows should have the specified width ({} columns)",
            WIDTH
        );

        // 3. Test grid content (all elements should be the default color code 0)
        let default_code = 0;
        let all_default = pg
            .grid_data
            .iter()
            .all(|row| row.iter().all(|&code| code == default_code));
        assert!(
            all_default,
            "All cells should be filled with the default color code {}",
            default_code
        );
    }

    #[test]
    /// Tests the dimensions method on a populated 10x10 grid.
    fn test_dimensions_valid_grid() {
        let pix_grid = get_test_pixgrid(); // 10 rows x 10 columns

        let (width, height) = pix_grid.dimensions();

        assert_eq!(width, 10, "Width should be 10 for the test grid");
        assert_eq!(height, 10, "Height should be 10 for the test grid");
    }

    #[test]
    /// Tests the dimensions method on an empty grid (no rows).
    fn test_dimensions_empty_grid() {
        // Create an empty PixGrid manually
        let pix_grid = PixGrid {
            color_map: HashMap::new(),
            grid_data: Vec::new(), // Empty grid_data
            cell_size: 10,
            grid_color: None,
        };

        let (width, height) = pix_grid.dimensions();

        assert_eq!(width, 0, "Width should be 0 for an empty grid");
        assert_eq!(height, 0, "Height should be 0 for an empty grid");
    }

    #[test]
    /// Tests the dimensions method on a grid that has rows but zero columns (should still return 0 width).
    fn test_dimensions_zero_width_rows() {
        // Create a PixGrid with 2 rows, each having 0 columns
        let pix_grid = PixGrid {
            color_map: HashMap::new(),
            grid_data: vec![Vec::new(), Vec::new()], // 2 rows, 0 columns each
            cell_size: 10,
            grid_color: None,
        };

        let (width, height) = pix_grid.dimensions();

        assert_eq!(height, 2, "Height should be 2 (number of rows)");
        assert_eq!(width, 0, "Width should be 0 (length of the first row)");
    }

    #[test]
    /// Tests the random_sparse_fill method to verify boundaries (0.0 and 1.0) and the distribution
    /// for an intermediate density.
    fn test_random_sparse_fill() {
        // 1. Setup: Create a 10x10 grid (100 total cells)
        let mut pg = PixGrid::new(10, 10);
        let total_cells = 100;
        const CELL_STATE: u8 = 1;

        // --- 1. Test density = 0.0 (Edge case: No activation) ---
        pg.random_sparse_fill(CELL_STATE, 0.0);
        let count_0 = pg
            .grid_data
            .iter()
            .flatten()
            .filter(|&&c| c == CELL_STATE)
            .count();
        assert_eq!(
            count_0, 0,
            "With a density of 0.0, no cells should be activated."
        );

        // --- 2. Test density = 1.0 (Edge case: All cells activated) ---
        pg.random_sparse_fill(CELL_STATE, 1.0);
        let count_1 = pg
            .grid_data
            .iter()
            .flatten()
            .filter(|&&c| c == CELL_STATE)
            .count();
        assert_eq!(
            count_1, total_cells,
            "With a density of 1.0, all cells should be activated."
        );

        // --- 3. Test intermediate density (0.35) ---
        const TARGET_DENSITY: f64 = 0.35;
        const TARGET_COUNT: usize = 35;

        // Define a tolerance to account for the random nature (e.g., +/- 15% of 100 cells)
        const TOLERANCE: usize = 15;

        // Fill the grid with the target density
        pg.random_sparse_fill(CELL_STATE, TARGET_DENSITY);

        let actual_count = pg
            .grid_data
            .iter()
            .flatten()
            .filter(|&&c| c == CELL_STATE)
            .count();

        // Define the boundaries of the expected range
        let lower_bound = TARGET_COUNT.saturating_sub(TOLERANCE); // saturating_sub ensures the lower bound is not negative
        let upper_bound = TARGET_COUNT + TOLERANCE;

        // A. Verify the range
        assert!(
            actual_count >= lower_bound && actual_count <= upper_bound,
            "The actual count ({}) should be in the range [{}, {}] for a density of {}",
            actual_count,
            lower_bound,
            upper_bound,
            TARGET_DENSITY
        );

        // B. Verify that the grid contains a mix of both codes (to ensure the random fill worked)
        let count_default = pg.grid_data.iter().flatten().filter(|&&c| c == 0).count();
        assert!(
            actual_count > 0,
            "At least one cell should be activated (code 1)."
        );
        assert!(
            count_default > 0,
            "At least one cell should be in the default state (code 0)."
        );
    }

    #[test]
    /// Tests the random grid filling method by checking proportion validation and the accuracy of final counts.
    fn test_random_fill_with_proportions() {
        // 1. Setup: Create a 10x10 grid (100 cells total)
        let mut pg = PixGrid::new(10, 10);
        let total_cells = 100;

        // Define target proportions:
        let valid_proportions = vec![
            (0, 0.10), // Target: 10 cells
            (1, 0.60), // Target: 60 cells
            (2, 0.30), // Target: 30 cells
        ];

        let invalid_proportions = vec![
            (0, 0.10),
            (1, 0.60),
            (2, 0.40), // Sums to 1.10 (Invalid)
        ];

        // --- 1. Failure Test: Invalid Proportions ---
        let result_invalid = pg.random_fill_with_proportions(&invalid_proportions);
        assert!(
            result_invalid.is_err(),
            "Should return an error if the sum of proportions is incorrect."
        );

        // --- 2. Success Test: Filling and Count Verification ---

        // Reset the grid and attempt valid fill
        let result_valid = pg.random_fill_with_proportions(&valid_proportions);
        assert!(
            result_valid.is_ok(),
            "Filling with valid proportions should succeed."
        );

        // A. Count the occurrences of each code after filling
        let mut counts: [usize; 3] = [0, 0, 0]; // For codes 0, 1, 2

        for row in &pg.grid_data {
            for &code in row {
                if (code as usize) < counts.len() {
                    counts[code as usize] += 1;
                }
            }
        }

        // B. Verify the total number of cells
        assert_eq!(
            counts.iter().sum::<usize>(),
            total_cells,
            "The total number of filled cells must equal 100."
        );

        // C. Verify individual counts against the target with tolerance

        // We expect the counts to be exactly the targets due to the `round()` and size adjustment logic,
        // but a small tolerance (e.g., 1) is used here to be robust against extremely minor floating-point quirks,
        // although generally unnecessary with a 100-cell grid.
        let tolerance = 1;

        // Target for Code 0 (10%)
        let target_0 = (total_cells as f64 * 0.10).round() as usize; // 10
        assert!(
            (counts[0] as isize - target_0 as isize).abs() <= tolerance as isize,
            "Code 0 count should be close to {} (Found: {})",
            target_0,
            counts[0]
        );

        // Target for Code 1 (60%)
        let target_1 = (total_cells as f64 * 0.60).round() as usize; // 60
        assert!(
            (counts[1] as isize - target_1 as isize).abs() <= tolerance as isize,
            "Code 1 count should be close to {} (Found: {})",
            target_1,
            counts[1]
        );

        // Target for Code 2 (30%)
        let target_2 = (total_cells as f64 * 0.30).round() as usize; // 30
        assert!(
            (counts[2] as isize - target_2 as isize).abs() <= tolerance as isize,
            "Code 2 count should be close to {} (Found: {})",
            target_2,
            counts[2]
        );
    }

    #[test]
    /// Tests the get_cell method for valid and invalid positions.
    fn test_get_cell() {
        let pix_grid = get_test_pixgrid(); // 10x10 grid

        // Définitions locales
        let center = CellPosition { x: 4, y: 4 };
        let valid_corner = CellPosition { x: 0, y: 0 };
        // Le test grid fait 10 de large/haut, donc 10 est hors limites (0..9)
        let out_of_bounds = CellPosition { x: 10, y: 10 };

        // 1. Valid position check (Center is code 0)
        assert_eq!(
            pix_grid.get_cell(center),
            Some(0),
            "Should retrieve code 0 for center cell"
        );

        // 2. Valid position check (Corner is code 0)
        assert_eq!(
            pix_grid.get_cell(valid_corner),
            Some(0),
            "Should retrieve code 0 for corner cell (0, 0)"
        );

        // 3. Invalid position check
        assert_eq!(
            pix_grid.get_cell(out_of_bounds),
            None,
            "Should return None for out-of-bounds position"
        );
    }

    #[test]
    /// Tests the set_cell method for valid and invalid positions, and verifies the change.
    fn test_set_cell() {
        let mut pix_grid = get_test_pixgrid(); // 10x10 grid, center is 0

        // Définitions locales
        let center = CellPosition { x: 4, y: 4 };
        let out_of_bounds = CellPosition { x: 100, y: 100 };
        const NEW_CODE: u8 = 99;

        // 1. Set cell at a valid position (CENTER)
        let result = pix_grid.set_cell(center, NEW_CODE);
        assert!(
            result.is_ok(),
            "Setting cell at a valid position should succeed"
        );

        // 2. Verify the change using get_cell
        assert_eq!(
            pix_grid.get_cell(center),
            Some(NEW_CODE),
            "get_cell should reflect the change to NEW_CODE"
        );

        // 3. Set cell at an invalid position (OUT_OF_BOUNDS)
        let error_result = pix_grid.set_cell(out_of_bounds, NEW_CODE);
        assert!(
            error_result.is_err(),
            "Setting cell out of bounds should return an error"
        );

        // 4. Verify that the grid data size hasn't changed (no panic/corruption)
        let (width, height) = pix_grid.dimensions();
        assert_eq!(width, 10);
        assert_eq!(height, 10);
    }

    #[test]
    /// Tests the get_moore_neighbors method for a central cell and a corner cell.
    fn test_get_moore_neighbors() {
        let mut pg = PixGrid::new(3, 3);

        // Define positions for testing
        let center = CellPosition { x: 1, y: 1 };
        let corner = CellPosition { x: 0, y: 0 };
        let out_of_bounds = CellPosition { x: 5, y: 5 };

        // Set codes for verification (creating a pattern of 1s around the center)
        // The grid will look like:
        // 0 1 0
        // 1 1 1
        // 0 1 0
        pg.set_cell(center, 1).unwrap();
        pg.set_cell(CellPosition { x: 0, y: 1 }, 1).unwrap(); // Left
        pg.set_cell(CellPosition { x: 2, y: 1 }, 1).unwrap(); // Right
        pg.set_cell(CellPosition { x: 1, y: 0 }, 1).unwrap(); // Top
        pg.set_cell(CellPosition { x: 1, y: 2 }, 1).unwrap(); // Bottom

        // --- 1. Test the center cell (should have 8 neighbors) ---
        let neighbors_center = pg.get_moore_neighbors(center);
        assert_eq!(
            neighbors_center.len(),
            8,
            "The center should have 8 neighbors"
        );

        // Count codes (4 zeros (diagonals), 4 ones (cardinals))
        let count_zeros = neighbors_center
            .iter()
            .filter(|&(_, code)| *code == 0)
            .count();
        let count_ones = neighbors_center
            .iter()
            .filter(|&(_, code)| *code == 1)
            .count();
        assert_eq!(
            count_zeros, 4,
            "The center should be surrounded by 4 zeros (diagonals)"
        );
        assert_eq!(
            count_ones, 4,
            "The center should be surrounded by 4 ones (cardinals)"
        );

        // --- 2. Test a corner cell (should have 3 neighbors) ---
        let neighbors_corner = pg.get_moore_neighbors(corner);
        assert_eq!(
            neighbors_corner.len(),
            3,
            "The corner should have 3 neighbors"
        );

        // Verify the positions of the corner's neighbors: (0, 1), (1, 0), (1, 1)
        let expected_positions = vec![
            CellPosition { x: 0, y: 1 }, // Bottom
            CellPosition { x: 1, y: 0 }, // Right
            CellPosition { x: 1, y: 1 }, // Diagonal
        ];
        for (pos, _) in neighbors_corner.iter() {
            assert!(
                expected_positions.contains(pos),
                "Neighbor position {:?} is not expected for the corner.",
                pos
            );
        }

        // --- 3. Test an out-of-bounds position (should be empty) ---
        let neighbors_oob = pg.get_moore_neighbors(out_of_bounds);
        assert!(
            neighbors_oob.is_empty(),
            "An out-of-bounds position should return no neighbors"
        );
    }

    #[test]
    /// Tests the get_von_neumann_neighbors method for a central cell and a corner cell.
    fn test_get_von_neumann_neighbors() {
        let mut pg = PixGrid::new(3, 3); // Creates a 3x3 grid (all codes are 0 by default)

        // Set codes for verification: Center and Left/Right are 1s, Top/Bottom are 2s.
        // The grid will look like:
        // 0 2 0
        // 1 1 1
        // 0 2 0
        let center = CellPosition { x: 1, y: 1 };
        let corner = CellPosition { x: 0, y: 0 };
        let out_of_bounds = CellPosition { x: 5, y: 5 };

        // Set center and cardinal neighbors
        pg.set_cell(center, 1).unwrap();
        pg.set_cell(CellPosition { x: 0, y: 1 }, 1).unwrap(); // Left (Code 1)
        pg.set_cell(CellPosition { x: 2, y: 1 }, 1).unwrap(); // Right (Code 1)
        pg.set_cell(CellPosition { x: 1, y: 0 }, 2).unwrap(); // Top (Code 2)
        pg.set_cell(CellPosition { x: 1, y: 2 }, 2).unwrap(); // Bottom (Code 2)

        // --- 1. Test the center cell (should have 4 neighbors) ---
        let neighbors_center = pg.get_von_neumann_neighbors(center);
        assert_eq!(
            neighbors_center.len(),
            4,
            "The center must have 4 cardinal neighbors"
        );

        // Check code counts: 2 cells with code 1 (Left/Right), 2 cells with code 2 (Top/Bottom)
        let count_ones = neighbors_center
            .iter()
            .filter(|&(_, code)| *code == 1)
            .count();
        let count_twos = neighbors_center
            .iter()
            .filter(|&(_, code)| *code == 2)
            .count();
        let count_zeros = neighbors_center
            .iter()
            .filter(|&(_, code)| *code == 0)
            .count();

        assert_eq!(count_ones, 2, "Should find 2 neighbors with code 1");
        assert_eq!(count_twos, 2, "Should find 2 neighbors with code 2");
        assert_eq!(count_zeros, 0, "Should find 0 neighbors with code 0");

        // --- 2. Test a corner cell (should have 2 neighbors) ---
        let neighbors_corner = pg.get_von_neumann_neighbors(corner);
        assert_eq!(
            neighbors_corner.len(),
            2,
            "The corner must have 2 cardinal neighbors"
        );

        // Expected neighbors for (0, 0) are (0, 1) [Bottom] and (1, 0) [Right]
        let expected_positions = vec![CellPosition { x: 0, y: 1 }, CellPosition { x: 1, y: 0 }];
        let found_positions: Vec<CellPosition> =
            neighbors_corner.iter().map(|(pos, _)| *pos).collect();

        assert!(found_positions.contains(&expected_positions[0]));
        assert!(found_positions.contains(&expected_positions[1]));

        // --- 3. Test an out-of-bounds position (should be empty) ---
        let neighbors_oob = pg.get_von_neumann_neighbors(out_of_bounds);
        assert!(
            neighbors_oob.is_empty(),
            "An out-of-bounds position should return no neighbors"
        );
    }
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
    fn test_export_pg_with_grid_color() {
        // Logic for creating the temporary path, moved from helper function
        let mut output_path = env::temp_dir();
        output_path.push("test_grid_with_color.pg");
        let path = output_path.as_path();

        // 1. Setup Data
        let mut color_map = HashMap::new();
        // Note: Adding colors in reverse order to test the sorting logic
        color_map.insert(2, Rgb([200, 100, 50]));
        color_map.insert(1, Rgb([255, 255, 255]));
        color_map.insert(0, Rgb([0, 0, 0]));

        let grid = PixGrid {
            color_map,
            grid_data: vec![vec![0, 0, 1, 1], vec![0, 2, 2, 1], vec![1, 2, 0, 0]],
            cell_size: 16,
            grid_color: Some(Rgb([50, 50, 50])),
        };

        // 2. Define Expected Output
        let expected_content = "\
cell_size = 16
grid_color = 50, 50, 50
0 = 0, 0, 0
1 = 255, 255, 255
2 = 200, 100, 50
---
0 0 1 1
0 2 2 1
1 2 0 0
";

        // 3. Call Function & Assert Success
        grid.export_pg(path)
            .expect("Failed to export PixGrid to .pg file");

        // 4. Read File Content and Assert (Logic moved from read_file_content)
        let mut file = File::open(path).expect("Failed to open exported file for reading");
        let mut actual_content = String::new();
        file.read_to_string(&mut actual_content)
            .expect("Failed to read exported file content to string");

        // Assert the content matches the expected string.
        assert_eq!(expected_content, actual_content);

        // 5. Cleanup
        fs::remove_file(path).expect("Failed to clean up temporary file");
    }

    #[test]
    fn test_export_pg_without_grid_color() {
        // Logic for creating the temporary path, moved from helper function
        let mut output_path = env::temp_dir();
        output_path.push("test_grid_without_color.pg");
        let path = output_path.as_path();

        // 1. Setup Data
        let mut color_map = HashMap::new();
        color_map.insert(5, Rgb([10, 20, 30]));

        let grid = PixGrid {
            color_map,
            grid_data: vec![vec![5, 5], vec![5, 5]],
            cell_size: 8,
            grid_color: None, // Test case: None
        };

        // 2. Define Expected Output
        let expected_content = "\
cell_size = 8
# grid_color = <R, G, B> # Omitted: No grid lines are drawn.
5 = 10, 20, 30
---
5 5
5 5
";

        // 3. Call Function & Assert Success
        grid.export_pg(path)
            .expect("Failed to export PixGrid without grid color");

        // 4. Read File Content and Assert (Logic moved from read_file_content)
        let mut file = File::open(path).expect("Failed to open exported file for reading");
        let mut actual_content = String::new();
        file.read_to_string(&mut actual_content)
            .expect("Failed to read exported file content to string");

        assert_eq!(expected_content, actual_content);

        // 5. Cleanup
        fs::remove_file(path).expect("Failed to clean up temporary file");
    }

    #[test]
    /// Tests that PNG generation creates a file with the correct dimensions.
    fn test_export_png_output() {
        let pix_grid = get_test_pixgrid();
        let output_path = PathBuf::from("temp_test_output.png");

        // Ensure cleanup of previous runs and defer cleanup for this run
        let _ = fs::remove_file(&output_path);

        // 1. Generate the PNG file
        pix_grid
            .export_png(&output_path)
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
            Err(e) => panic!("Failed to open exported PNG for verification: {}", e),
        }

        // 3. Cleanup
        fs::remove_file(&output_path).expect("Failed to clean up exported PNG file");
    }

    #[test]
    /// Tests that SVG generation creates a file with the correct header and stroke color.
    fn test_export_svg_output() {
        let pix_grid = get_test_pixgrid();
        let output_path = PathBuf::from("temp_test_output.svg");

        // Ensure cleanup of previous runs and defer cleanup for this run
        let _ = fs::remove_file(&output_path);

        // 1. Generate the SVG file
        pix_grid
            .export_svg(&output_path)
            .expect("SVG generation failed");

        // 2. Assert file was created
        assert!(
            output_path.exists(),
            "SVG file should exist after generation"
        );

        // 3. Assert content (400x400) and the stroke color
        let content = fs::read_to_string(&output_path).expect("Failed to read exported SVG file");

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
        fs::remove_file(&output_path).expect("Failed to clean up exported SVG file");
    }
}
