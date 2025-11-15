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
