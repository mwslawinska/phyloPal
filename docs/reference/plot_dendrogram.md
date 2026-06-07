# Plot a dendrogram as ggplot with optional tip aesthetics (color/fill/shape/size) and labels

Converts a `dendrogram` into ggplot layers using `ggdendro`. Draws
branches as line segments and overlays tip (leaf) points. Tip aesthetics
can be mapped to outline color (`color`), interior fill (`fill`; visible
for shapes 21-25), shape (`shape`), and size (`size`).

## Usage

``` r
plot_dendrogram(
  dend,
  metadata = NULL,
  label_from = NULL,
  color_by = NULL,
  color_palette = NULL,
  fill_by = NULL,
  fill_palette = NULL,
  shape_by = NULL,
  shape_palette = NULL,
  size_by = NULL,
  size_palette = NULL,
  tip_meta = NULL,
  label_col = "label",
  point_size = 2,
  point_stroke = NULL,
  default_shape = 21,
  orientation = c("top", "right", "bottom", "left"),
  first_branch = c("bottom", "top"),
  add_tip_labels = FALSE,
  tip_label_size = 3,
  leaf_order = NULL,
  rotate_to_leaf_order = TRUE,
  theme_obj = ggplot2::theme_void(),
  guides_obj = NULL,
  legend_position = "right",
  legend_position_inside = NULL
)
```

## Arguments

- dend:

  A `dendrogram` object, or a list from
  [`build_dendrogram`](build_dendrogram.md) containing `$dendrogram`.

- metadata:

  Optional data frame used to build `tip_meta` automatically
  (recommended).

- label_from:

  Column in `metadata` that matches the dendrogram leaf labels.

- color_by:

  Optional column mapped to point outline color.

- color_palette:

  Optional named vector for `color_by` levels.

- fill_by:

  Optional column mapped to point interior fill.

- fill_palette:

  Optional named vector for `fill_by` levels.

- shape_by:

  Optional column mapped to point shape.

- shape_palette:

  Optional values passed to
  [`ggplot2::scale_shape_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

- size_by:

  Optional column mapped to point size.

- size_palette:

  Optional named numeric vector passed to
  [`ggplot2::scale_size_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

- tip_meta:

  Optional precomputed tip metadata. Must contain `label_col` and any
  requested mapping columns.

- label_col:

  Name of the label column inside `tip_meta` (default: `"label"`).

- point_size:

  Numeric tip point size (default: 2). Ignored when `size_by` is used.

- point_stroke:

  Numeric outline width for points (`stroke`); only applies to shapes
  21–25. If `NULL`, uses ggplot2 default.

- default_shape:

  Integer shape used when `fill_by` is set but `shape_by` is `NULL`
  (default: 21). This ensures mapped fills are visible.

- orientation:

  Dendrogram orientation:

  `"top"`

  :   Tips at bottom (default ggdendro layout).

  `"bottom"`

  :   Tips at top (vertical mirror).

  `"right"`

  :   Tips on left, branches extend to the right (rotated 90°).

  `"left"`

  :   Tips on right, branches extend to the left (rotated 90° +
      mirrored).

- first_branch:

  For `orientation="left"` or `"right"`, controls whether the first leaf
  in the layout appears at the top (`"top"`) or bottom (`"bottom"`).
  This is implemented by reversing the (post-flip) vertical ordering via
  [`scale_x_reverse()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html).
  For `orientation="top"` or `"bottom"`, this parameter is currently
  ignored.

- add_tip_labels:

  Logical; if `TRUE`, draw tip labels.

- tip_label_size:

  Numeric text size for tip labels.

- leaf_order:

  Optional character vector giving the desired leaf order (matching
  dendrogram labels).

- rotate_to_leaf_order:

  Logical; if `TRUE` (default), rotate dendrogram to match `leaf_order`
  using
  [`dendextend::rotate()`](https://talgalili.github.io/dendextend/reference/rotate.html)
  (requires dendextend).

- theme_obj:

  ggplot2 theme to apply (default:
  [`ggplot2::theme_void()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).

- guides_obj:

  Optional `ggplot2::guides(...)` added to the plot.

- legend_position:

  Legend position passed to `theme(legend.position=...)`.

- legend_position_inside:

  Numeric length-2 vector (x, y) between 0 and 1 used when
  `legend_position="inside"`.

## Value

A `ggplot` object. Leaf order (left-to-right) is stored as
`attr(plot, "leaf_order")`.

## Details

Leaf order (left-to-right in the rendered dendrogram layout *before
orientation transforms*) is stored on the returned plot as
`attr(plot, "leaf_order")`. This is useful to order the alluvial axis to
match the dendrogram.
