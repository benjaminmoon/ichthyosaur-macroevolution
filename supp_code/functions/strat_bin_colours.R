#!/usr/bin/env Rscript

# This script includes a function to more quickly input RGB colours in a 0–255
# scale, and a list of RGB colour values for the Mesozoic; periods, epochs, and
# ages.
# From Moon & Stubbs <Title>, <doi>

# Contents:
#   - RBG
#   - bin_colours

# function to always use 256 steps (0–255)
RGB <- function(red, green, blue) {
  # Converts the given red, green, and blue values into a RGB hex code.
  #
  # Args:
  #   red, green, blue: each of red, green, and blue values on a 0–255 scale.
  #
  # Returns: a RGB hex code for e.g. plotting with.
  rgb(red, green, blue, maxColorValue = 255)
}

# create Mesozoic colours list
bin_colours <- list(periods = list(Triassic   = RGB(129,  43, 153),
                                   Jurassic   = RGB( 52, 178, 201),
                                   Cretaceous = RGB(127, 198,  78)),
                    epochs  = list(Early_Triassic   = RGB(152,  57, 153),
                                   Middle_Triassic  = RGB(177, 104, 177),
                                   Late_Triassic    = RGB(189, 140, 195),
                                   Early_Jurassic   = RGB( 66, 174, 208),
                                   Middle_Jurassic  = RGB(128, 207, 216),
                                   Late_Jurassic    = RGB(179, 227, 238),
                                   Early_Cretaceous = RGB(140, 205,  87),
                                   Late_Cretaceous  = RGB(166, 216,  74)),
                    ages    = list(Induan        = RGB(164,  70, 159),
                                   Olenekian     = RGB(176,  81, 165),
                                   Anisian       = RGB(188, 117, 183),
                                   Ladinian      = RGB(201, 131, 191),
                                   Carnian       = RGB(201, 155, 203),
                                   Norian        = RGB(214, 170, 211),
                                   Rhaetian      = RGB(227, 185, 219),
                                   Hettangian    = RGB( 78, 179, 211),
                                   Sinemurian    = RGB(103, 188, 216),
                                   Pliensbachian = RGB(128, 197, 221),
                                   Toarcian      = RGB(153, 206, 227),
                                   Aalenian      = RGB(154, 217, 221),
                                   Bajocian      = RGB(166, 221, 224),
                                   Bathonian     = RGB(179, 226, 227),
                                   Callovian     = RGB(191, 231, 229),
                                   Oxfordian     = RGB(191, 231, 241),
                                   Kimmeridgian  = RGB(204, 236, 244),
                                   Tithonian     = RGB(217, 241, 247),
                                   Berriasian    = RGB(140, 205,  96),
                                   Valanginian   = RGB(153, 211, 106),
                                   Hauterivian   = RGB(166, 217, 117),
                                   Barremian     = RGB(179, 223, 127),
                                   Aptian        = RGB(191, 228, 138),
                                   Albian        = RGB(204, 234, 151),
                                   Cenomanian    = RGB(179, 222,  83),
                                   Turonian      = RGB(191, 227,  93),
                                   Coniacian     = RGB(204, 233, 104),
                                   Santonian     = RGB(217, 239, 116),
                                   Campanian     = RGB(230, 244, 127),
                                   Maastrichtian = RGB(242, 250, 140)))