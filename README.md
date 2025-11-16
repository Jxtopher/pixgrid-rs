# PixGrid

The PixGrid library is a utility written in Rust designed to convert simple, human-readable text definitions of a pixel grid into high-quality image formats (PNG or SVG).

This format is ideal for defining pixel art, simple maps, or custom icons in a clear, text-based manner.


## PixGrid Format Specification (.pg)

The PixGrid file format (typically using the `.pg` extension) is divided into three main sections. The parser ignores lines starting with `#` (comments).

### 1. Rendering Parameters

This section defines the structural properties of the output image.

| Format | Description |
|--------|-------------|
| cell_size | The size (width and height) in pixels of each individual grid cell in the final image. Default is 10.|
| grid_color | Defines the color of the 1-pixel wide grid lines drawn between cells. If this parameter is omitted, no grid lines are drawn. |

### 2. Color Definitions

| Format | Description |
|--------|-------------|
| `[Code] = R, G, B` | Defines the color map entry. `Code` must be a single `u8` integer. `R`, `G`, `B` must be comma-separated `u8` (0-255) values. |

### 3. Grid Data

This section contains the actual pixel layout, separated from the configuration by ---.

|   Format  | Description |
|-----------|-------------|
| `---`     |Mandatory separator marking the start of the grid data.|
| `[codes]` | Each line represents a row in the grid. Color codes must be separated by whitespace and must refer to a defined entry in the color map. |


## Example PixGrid File

This example file defines a simple 10x10 pattern with cells of 40x40 pixels, using black lines to delineate the grid.

```txt
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
```

## Snippet

This small Rust program is a command-line tool that uses the `pixgrid` library to visualize data.

1. It reads a plain text definition file (`.pg`) that contains color codes and grid layout information.
2. It uses `PixGrid::parse` to convert this text into a structured grid object.
3. Based on the output file's extension (`.png` or `.svg`), it calls the corresponding method (`generate_png` or `generate_svg`) to draw the pixel art image and save it to the disk.


```rust
use pixgrid::PixGrid;
use std::error::Error;
use std::fs;
use std::path::Path;
fn main() -> Result<(), Box<dyn Error>> {
    let input_file_path = "instances/nested_squares.pg";
    // Modify here to choose the output extension and format.
    let output_file_name = "nested_squares.svg";
    let output_path = Path::new(output_file_name);
    // 1. Reading the file
    println!("Reading file: {}", input_file_path);
    let input_data = fs::read_to_string(input_file_path).map_err(|e| {
        format!(
            "Error reading file {}: {}",
            input_file_path, e
        )
    })?;
    // 2. Parsing colors and the grid (PixGrid contains everything)
    // Using the static method PixGrid::parse
    let pg = PixGrid::parse(&input_data)?;
    println!(
        "Defined colors: {:?}",
        pg.color_map.keys().collect::<Vec<_>>()
    );
    println!(
        "Grid read: {} rows x {} columns",
        pg.grid_data.len(),
        pg.grid_data.first().map_or(0, |r| r.len())
    );
    // 3. Image generation (selection logic)
    let extension = output_path
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .unwrap_or_default();
    match extension {
        "png" => {
            println!("Generating in PNG format...");
            // Using the instance method
            pg.generate_png(output_path)?;
        }
        "svg" => {
            println!("Generating in SVG format...");
            // Using the instance method
            pg.generate_svg(output_path)?;
        }
        _ => {
            return Err(format!(
                "Unsupported file extension: '{}'. Use '.png' or '.svg'.",
                extension
            )
            .into());
        }
    }
    println!("Generation finished successfully!");
    Ok(())
}
```
