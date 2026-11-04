# Analysis design schematic: GBC GWAS hit -> mediator -> gallbladder cancer,
# with the path labels used in docs/methods.md.
#
# Usage: Rscript scripts/design-figure.R
# Writes figures/analysis-design.pdf and figures/analysis-design.png.

library(grid)

draw_design <- function() {
  grid.newpage()
  box_w <- unit(1.25, "inches")
  box_h <- unit(0.8, "inches")
  g_size <- unit(0.45, "inches")
  arrow_head <- arrow(length = unit(0.12, "inches"), type = "closed")
  line_gp <- gpar(lwd = 2, col = "black", fill = "black")

  node <- function(x, y, label, w = box_w, h = box_h, fontsize = 13) {
    grid.roundrect(x, y, w, h, r = unit(0.08, "snpc"), gp = gpar(lwd = 2, fill = "white"))
    grid.text(label, x, y, gp = gpar(fontsize = fontsize, fontface = "bold"))
  }
  path_label <- function(label, x, y) {
    grid.text(label, x, y, gp = gpar(fontsize = 16, fontface = "bold.italic"))
  }

  hit <- c(0.15, 0.55)
  med <- c(0.5, 0.55)
  gbc <- c(0.85, 0.55)
  half_w <- 0.5 * convertWidth(box_w, "npc", valueOnly = TRUE)
  half_h <- 0.5 * convertHeight(box_h, "npc", valueOnly = TRUE)

  # d: residual direct path from the GBC hit to GBC that bypasses the mediator.
  grid.curve(hit[1], med[2] + half_h, gbc[1], med[2] + half_h, curvature = -0.45,
             ncp = 20, square = FALSE, arrow = arrow_head, gp = line_gp)
  path_label("d", 0.5, 0.94)

  # a: GBC hit -> mediator; b: mediator -> GBC.
  grid.segments(hit[1] + half_w, hit[2], med[1] - half_w - 0.005, med[2],
                arrow = arrow_head, gp = line_gp)
  path_label("a", mean(c(hit[1], med[1])), hit[2] + 0.05)
  grid.segments(med[1] + half_w, med[2], gbc[1] - half_w - 0.005, gbc[2],
                arrow = arrow_head, gp = line_gp)
  path_label("b", mean(c(med[1], gbc[1])), med[2] + 0.05)

  # Other mediator instruments, used to estimate b.
  g_x <- c(0.35, 0.45, 0.55, 0.65)
  g_y <- 0.17
  g_half <- 0.5 * convertHeight(g_size, "npc", valueOnly = TRUE)
  for (x in g_x) {
    grid.segments(x, g_y + g_half, med[1] + (x - med[1]) * 0.25, med[2] - half_h - 0.01,
                  arrow = arrow_head, gp = line_gp)
    node(x, g_y, "G", g_size, g_size, fontsize = 12)
  }
  grid.text("Other mediator instruments", 0.5, 0.05, gp = gpar(fontsize = 13, fontface = "bold"))

  node(hit[1], hit[2], "GBC\nGWAS hit")
  node(med[1], med[2], "Mediator")
  node(gbc[1], gbc[2], "Gallbladder\ncancer")
}

dir.create("figures", showWarnings = FALSE)
pdf("figures/analysis-design.pdf", width = 8, height = 5)
draw_design()
invisible(dev.off())
png("figures/analysis-design.png", width = 8, height = 5, units = "in", res = 300)
draw_design()
invisible(dev.off())
message("Wrote figures/analysis-design.pdf and figures/analysis-design.png")
