
#' View the colors of a HTML list
#'
#' @param html_codes HTML codes you want to have a look at
#'
#' @return A figure showing the palette colors
#' @export
#'
#' @examples
#' colors_vec <- ggsci::pal_futurama()(4)
#'
#' show_palette_colors(colors_vec)
#'
show_palette_colors <- function(html_codes){
  image(1:length(html_codes), 1, as.matrix(1:length(html_codes)),
        col = html_codes,
        xlab = "", ylab = "")}
