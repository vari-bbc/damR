
#' Produce volcano plot data for a single contrast
#'
#' @inheritParams PlotVolcano
#'
#' @returns data.table ready for plotting with `VolcanoPlotter()`
#'
#' @importFrom collapse fsubset fmutate
#' @importFrom data.table fcase
#'
#' @export
VolcanoData <- function(vpd, contrastGroups, contrastTitle, alphaLevel = 0.05, cap = 10) {
  contrastStr <- MakeTheContrasts(contrastGroups, contrastTitle)
  fsubset(vpd, Contrast == contrastStr) |>
    fmutate(
      logFC = pmax(pmin(logFC, cap), -cap), # cap between -cap and cap
      Direction = fcase(
        logFC > 0 & adj.P.Val < alphaLevel , "Up regulated"   ,
        logFC < 0 & adj.P.Val < alphaLevel , "Down regulated" ,
        default = "NS"
      )
    )
}

#' Volcano plot
#'
#' @param theData The cleaned? fit? data
#' @param theAlpha Alpha level
#' @returns Volcano ggplot
#'
#' @importFrom ggiraph geom_point_interactive girafe girafe_options opts_selection
#' @importFrom ggplot2 aes facet_wrap geom_hline ggplot ggtitle labs scale_color_manual theme theme_minimal ylab
#' @importFrom rlang .data
#'
#' @export
MakeVolcano <- function(theData, theAlpha = 0.05) {
  maxAbsLogFC <- max(abs(theData$logFC), na.rm = TRUE)
  if (!is.finite(maxAbsLogFC) || maxAbsLogFC == 0) {
    maxAbsLogFC <- 1
  }

  ## ERW
  #dat_sort <- sort(theData$logFC)
  #gaps <- which(diff(dat_sort) == max(diff(dat_sort)))
  # if(dat_sort[gaps] < 0) { ## create break when the largest diff is positive or negative
  #   gaps_total <- dat_sort[gaps - -1:0]
  #
  # } else {
  #   gaps_total <- dat_sort[gaps + 0:1]
  # }

  # order data frame by logFC
  # dat_sort <- theData[order(theData$logFC), ]

  ## if the difference between any two consecutive points is > value, flag as outlier?
  # outliers <- dat_sort[diff(dat_sort$logFC) > 1000, ]

  ## then make plot:

  theData$mod_logFC <- theData$logFC
  theData$mod_logFC <- ifelse(theData$mod_logFC < -10, -10, theData$mod_logFC)
  theData$mod_logFC <- ifelse(theData$mod_logFC > 10, 10, theData$mod_logFC)

  volc <- ggplot(
    theData,
    aes(
      x = .data$mod_logFC,
      y = -log10(.data$adj.P.Val),
      color = .data$Direction,
      data_id = .data$Compound,
      tooltip = .data$Compound,
      label = .data$Compound
    )
  ) +
    geom_point_interactive(size = 2) +
    theme_minimal(7) +
    scale_color_manual(
      values = c(
        "Up regulated" = "#276419",
        "Down regulated" = "#8E0152",
        "NS" = "grey"
      ),
      drop = F
    ) +
    geom_hline(
      yintercept = -log10(theAlpha),
      linetype = "dotted",
      color = "red",
      linewidth = 0.5
    ) +
    labs(
      subtitle = theData$Contrast[1],
      caption = "all logFC capped at 10/-10 for plotting purposes"
    ) +
    # geom_text_interactive(alpha=0,size=4) +
    # geom_text_repel(theData = theData[theData$adj.P.Val < 0.05,],
    #                 aes(x = logFC,y = -log10(adj.P.Val), label = label),
    #                 max.overlaps = Inf) +
    ylab("-log10(FDR p)") +
    #coord_cartesian(xlim = c(-maxAbsLogFC, maxAbsLogFC)) + # ERW
    theme(legend.position = "none")

  # ## ERW if max diff is larger than 1000
  # if(max(diff(dat_sort)) > 10000) {
  #   volc <- volc + scale_x_break(breaks = gaps_total)
  # }

  # If multiple sub contrasts, facet wrap by them so they stay together
  if (length(unique(theData$Contrast)) > 1) {
    volc <- volc + facet_wrap(~Contrast)
  }

  return(volc)
}

#' Volcano plots for all contrasts
#'
#' @param theContrasts The contrasts for comparison
#' @param ebayesFit The fit data, padjusted?
#' @param theAlpha Alpha level
#' @returns A list of plots and data
#'
#' @importFrom ggiraph girafe_css opts_selection opts_sizing
#'
#' @export
MakeAllVolcanoes <- function(theContrasts, ebayesFit, theAlpha) {
  allVolcanoes <- list()
  for (i in 1:length(unlist(theContrasts))) {
    # Get data for volcano plot
    tops.all <- topTable(
      ebayesFit,
      adjust.method = "BH",
      p.value = 1,
      number = Inf,
      coef = i
    )

    # Add current contrast
    theData <- data.frame(Contrast = unlist(theContrasts)[i], tops.all)

    # Make sure rownames match and set Compound to be the rownames
    rownames(theData) <- rownames(tops.all)
    theData$Compound <- rownames(theData)

    # Classify significance
    theData$Significant <- ifelse(
      theData$adj.P.Val < theAlpha,
      "Significant",
      "NS"
    )

    # Add labels for the top 25
    theData$label <- ""
    theData <- theData[order(theData$adj.P.Val, decreasing = F), ]
    theData$label[1:25] <- theData$Compound[1:25]

    # Add direction for significant metabolites
    theData$Direction <- ifelse(
      theData$logFC > 0,
      "Up regulated",
      "Down regulated"
    )
    theData$Direction[theData$Significant == "NS"] <- "NS"
    theData$Direction <- factor(
      theData$Direction,
      levels = c("Down regulated", "NS", "Up regulated")
    )

    # Make volcanoes
    allVolcanoes[[i]] <- MakeVolcano(theData, theAlpha)
  }

  p <- list()

  lng <- 1:length(theContrasts[[1]])

  # Add titles to the plots based on title given at beginning
  # wrappedPlots <- wrap_plots(lapply(lng,function(x) allVolcanoes[[x]] + ggtitle(names(theContrasts)[1]) + theme(legend.position = 'top') ))
  # p[[1]] <- wrappedPlots
  p <- lapply(lng, function(x) {
    allVolcanoes[[x]] +
      ggtitle(names(theContrasts)[1]) +
      theme(legend.position = 'top')
  })

  # Combine all plots of all contrasts
  # g <- wrap_plots(p) + plot_layout(ncol = 1)
  g <- Reduce(`/`, p)

  # Turn ggplot into girafe
  x <- girafe(
    ggobj = g,
    height_svg = 4 * length(unlist(theContrasts)) #, ## ERW x <- girafe(ggobj = g) error in grid.newpage(): svgd only supports one page
    # options = list(
    #   opts_sizing(rescale = FALSE) # Stops it from filling the container
    # )
  )

  VolcanoPlots <- girafe_options(
    x,
    opts_selection(
      type = "single",
      only_shiny = FALSE,
      girafe_css(
        css = "fill:orange; stroke:black; 
                                                   r:3pt; font-style:italic; 
                                                   font-size:5px; max-height: 1000%;"
      )
    ),
    opts_sizing(rescale = TRUE)
  )

  return(list(VolcanoPlots, theData))
}
