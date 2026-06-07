# Add alpha transparency to hex colors

Add alpha transparency to hex colors

## Usage

``` r
add_alpha(hex, alpha = 1)
```

## Arguments

- hex:

  Character vector of hex colors

- alpha:

  Numeric value between 0 (transparent) and 1 (opaque)

## Value

Character vector of hex colors with alpha channel

## Examples

``` r
add_alpha("#FF0000", 0.5)  # Semi-transparent red
#> [1] "#FF00007F"
add_alpha(c("#FF0000", "#00FF00"), 0.7)
#> [1] "#FF0000B2" "#00FF00B2"
```
