# Base graphics exports: every pair is shown, with independent comparison scales.
format_probability <- function(x) ifelse(is.finite(x), formatC(x, digits = 2, format = "g"), "NA")
wrap_label <- function(x, width = 88) paste(strwrap(x, width), collapse = "\n")

draw_estimates <- function(b, lo, hi, labels, colours, title, subtitle, xlab) {
  available <- is.finite(b) & is.finite(lo) & is.finite(hi)
  limits <- range(c(0, lo[available], hi[available]))
  if (diff(limits) == 0) limits <- c(-1, 1)
  limits <- limits + c(-1, 1) * diff(limits) * 0.08
  y <- rev(seq_along(b))
  par(mar = c(2.7, 7, 3.8, 1), mgp = c(1.6, 0.45, 0), tcl = -0.2)
  plot(NA, xlim = limits, ylim = c(0.5, length(b) + 0.5), yaxt = "n",
    ylab = "", xlab = xlab, cex.lab = 0.75, cex.axis = 0.8, bty = "n")
  abline(v = 0, lty = 3, col = "grey65")
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.8)
  segments(lo[available], y[available], hi[available], y[available], col = colours[available], lwd = 1.7)
  points(b[available], y[available], col = colours[available], pch = 19, cex = 0.9)
  if (any(!available)) text(mean(limits), y[!available], "Unavailable (see CSV reason)",
    cex = 0.75, col = "grey40")
  mtext(wrap_label(title), side = 3, line = 1.7, adj = 0, cex = 0.85, font = 2)
  mtext(subtitle, side = 3, line = 0.1, adj = 0, cex = 0.7)
}

export_pages <- function(path, n, per_page, draw_page, width = 10, height = 11) {
  pages <- split(seq_len(n), ceiling(seq_len(n) / per_page))
  # Embed fonts so publication PDFs preserve spacing across viewers/platforms.
  grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  tryCatch(for (page in seq_along(pages)) draw_page(pages[[page]], page, length(pages)),
    finally = grDevices::dev.off())
  pngs <- character(length(pages))
  for (page in seq_along(pages)) {
    suffix <- if (page == 1) ".png" else sprintf("-page-%02d.png", page)
    pngs[page] <- sub("\\.pdf$", suffix, path)
    grDevices::png(pngs[page], width = width, height = height, units = "in", res = 150, type = "cairo")
    tryCatch(draw_page(pages[[page]], page, length(pages)), finally = grDevices::dev.off())
  }
  pngs
}

export_comparison_plot <- function(x, path) {
  export_pages(path, nrow(x), 6, function(rows, page, pages) {
    par(mfrow = c(6, 1), oma = c(2, 0, 3, 0))
    for (row in rows) {
      z <- x[row, ]
      title <- paste(z$mediator, "|", z$mediator_id, "|", z$target_snp)
      subtitle <- paste0("Effect allele: ", z$effect_allele, "; Q BH p=", format_probability(z$q_padj),
        "; ", ifelse(is.na(z$q_padj), "unavailable", ifelse(z$q_padj < 0.05, "incompatible", "unresolved")))
      draw_estimates(c(z$target_b, z$comparator_b), c(z$target_lo, z$comparator_lo),
        c(z$target_hi, z$comparator_hi),
        c(paste0("Exact hit (n=", z$target_nsnp, ")"), paste0("Other SNPs (n=", z$comparator_nsnp, ")")),
        c("#2166ac", "#b35806"), title, subtitle, "GBC log odds / saved mediator unit (panel-specific scale)")
    }
    mtext(sprintf("Exact hit versus other instruments | page %d/%d", page, pages), outer = TRUE,
      side = 3, line = 1, font = 2, cex = 1.1)
    mtext("95% CIs; primary +/-1 Mb exclusions; conservative palindrome removal. Q family: unique pairs.",
      outer = TRUE, side = 1, line = 0.5, cex = 0.75)
  })
}

export_decomposition_plot <- function(x, path) {
  colours <- c("#2166ac", "#1b7837", "#b35806")
  export_pages(path, nrow(x), 11, function(rows, page, pages) {
    z <- x[rows, ]
    b <- c(z$t, z$i, z$d); lo <- c(z$t_lo, z$i_lo, z$d_lo); hi <- c(z$t_hi, z$i_hi, z$d_hi)
    limits <- range(c(0, lo[is.finite(lo)], hi[is.finite(hi)]))
    if (diff(limits) == 0) limits <- c(-1, 1)
    limits <- limits + c(-1, 1) * diff(limits) * 0.07
    par(mfrow = c(1, 1), mar = c(5, 19, 5, 2), mgp = c(2.6, 0.7, 0))
    y <- rev(seq_len(nrow(z)))
    plot(NA, xlim = limits, ylim = c(0.4, nrow(z) + 0.9), yaxt = "n", ylab = "",
      xlab = "GBC log odds per recorded target effect allele", bty = "n")
    abline(v = 0, lty = 3, col = "grey60")
    labels <- vapply(seq_len(nrow(z)), function(k) {
      paste(paste(strwrap(z$mediator[k], 34), collapse = "\n"),
        paste0(z$target_snp[k], " | EA ", z$effect_allele[k]),
        paste0("Other n=", z$comparator_nsnp[k], " | residual BH p=", format_probability(z$d_padj[k])), sep = "\n")
    }, character(1))
    axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.68)
    for (component in seq_along(colours)) {
      field <- c("t", "i", "d")[component]
      available <- is.finite(z[[field]]) & is.finite(z[[paste0(field, "_lo")]]) & is.finite(z[[paste0(field, "_hi")]])
      yp <- y + c(0.22, 0, -0.22)[component]
      segments(z[[paste0(field, "_lo")]][available], yp[available],
        z[[paste0(field, "_hi")]][available], yp[available], col = colours[component], lwd = 1.5)
      points(z[[field]][available], yp[available], pch = 19, col = colours[component], cex = 0.75)
      if (any(!available)) text(limits[2], yp[!available], paste(c("Total", "Indirect", "Residual")[component], "unavailable"),
        adj = 1, cex = 0.6, col = colours[component])
    }
    title(main = sprintf("Model-based components | page %d/%d", page, pages), line = 3.2)
    legend("top", inset = c(0, -0.035), xpd = TRUE, horiz = TRUE, bty = "n", cex = 0.8,
      legend = c("Total", "Indirect", "Residual direct"), col = colours, pch = 19)
    mtext("95% delta-method CIs; zero-covariance working assumption; each mediator analysed separately.",
      side = 1, line = 3.6, cex = 0.7)
  }, width = 12)
}

export_pooled_plot <- function(x, path) {
  export_pages(path, nrow(x), 6, function(rows, page, pages) {
    par(mfrow = c(6, 1), oma = c(2, 0, 3, 0))
    for (row in rows) {
      z <- x[row, ]
      draw_estimates(z$b, z$lo, z$hi, paste0("Pooled (n=", z$nsnp, ")"), "#1b7837",
        paste(z$mediator, "|", z$mediator_id),
        paste0(z$method, "; BH p=", format_probability(z$pval_bh), "; excluded target/region n=",
          z$n_excluded_target + z$n_excluded_region),
        "GBC log odds / saved mediator unit (panel-specific scale)")
    }
    mtext(sprintf("Pooled mediator MR | page %d/%d", page, pages), outer = TRUE,
      side = 3, line = 1, font = 2, cex = 1.1)
    mtext("95% CIs; target-locus union excluded before harmonisation; BH family: pooled mediators.",
      outer = TRUE, side = 1, line = 0.5, cex = 0.75)
  })
}
