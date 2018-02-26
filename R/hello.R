# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, worlds!")
}


#' Add together two numbers.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' add(10, 1)
add <- function(x, y) {
    x + y
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Welcome to maser v 0.0.1")
}

my_fun <- function(a, b) {
    if (!requireNamespace("pkg", quietly = TRUE)) {
        stop("Package \"pkg\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
}

