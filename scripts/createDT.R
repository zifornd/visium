#' Create DT table with breaks.
#'
#' @param sample_name Name of sample.
#' @param res List of results tables.
#' @return list of header, dt table, and breaks.
createDT <- function(sample_name, res) {

  header <- htmltools::tags$h4(sample_name)

  dt <- DT::datatable(res[[sample_name]],
                      class = 'cell-border stripe',
                      extensions = 'Buttons',
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
                      filter = "top",
                      width = "100%")

  br1 <- htmltools::tags$br()
  br2 <- htmltools::tags$br()
  br3 <- htmltools::tags$br()
  br4 <- htmltools::tags$br()
  br5 <- htmltools::tags$br()
  br6 <- htmltools::tags$br()
  br7 <- htmltools::tags$br()
  br8 <- htmltools::tags$br()


  return(list(header, dt, br1, br2, br3, br4, br5, br6, br7, br8))
}

#' Create DT table with breaks.
#'
#' @param sample_name Name of sample.
#' @param res List of results tables.
#' @param unlist If we want to remove sample information and collapse into sample/cluster named list. 
#' @param html If we want to return html tag output for tagList rendering.
#' @param targets targets details for DT table JS.
#' @return DT datatable OR list of header, DT table, and breaks.
createDTFea <- function(sample_name, res, group,
                        unlist = FALSE, html = FALSE, targets = c(11, 12)) {

  # If we want to remove sample information and collapse into sample/cluster named list
  if (unlist) {

    res <- unlist(res, recursive = FALSE)

    # remove empty elements for clusters without any significant results
    resbind <- res[sapply(res, nrow) > 0]

    # bind into single results table
    resbind <- dplyr::bind_rows(resbind, .id = group)

  } else {

    # remove empty elements for clusters without any significant results
    resbind <- res[sapply(res[[sample_name]], nrow) > 0]

    # bind into single results table
    resbind <- dplyr::bind_rows(resbind, .id = group)
  }

  # get dt table object
  dt <- plotFeaDT(resbind, pageLength = 5, targets = c(11, 12))

  # get html tags from html tools if required for multiple dt tables.
  if (html) {

    header <- htmltools::tags$h4(sample_name)

    br1 <- htmltools::tags$br()
    br2 <- htmltools::tags$br()
    br3 <- htmltools::tags$br()
    br4 <- htmltools::tags$br()
    br5 <- htmltools::tags$br()
    br6 <- htmltools::tags$br()
    br7 <- htmltools::tags$br()
    br8 <- htmltools::tags$br()

    dt_complete <- list(header, dt, br1, br2, br3, br4, br5, br6, br7, br8)

  } else {

    dt_complete <- dt

  }

  return(dt_complete)

}

#' Create DT table with JS length restrictions for FEA.
#'
#' @param sample_name Results table following FEA.
#' @param pageLength pageLength.
#' @param targets targets for JS to shorten output displayed (suffix ...).
#' @return DT datatable.
plotFeaDT <- function(resbind, pageLength = 5, targets = c(11, 12)) {

    DT::datatable(resbind, class = 'cell-border stripe',
                  extensions = 'Buttons', rownames = FALSE,
                  options = list(
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                    pageLength = pageLength,
                    columnDefs = list(
                      list(targets = targets,
                           render = JS(
                             "function(data, type, row, meta) {",
                             "return type === 'display' && data.length > 6 ?",
                             "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
                             "}")))),
                  filter = "top")
}