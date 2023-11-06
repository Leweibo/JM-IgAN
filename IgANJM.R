
dynamicplotIgAN <- function(JMmodel = jointFit_clinical5 ,newD = ND,years_pred = 10,length.out = 50) {
  tmax = max(newD$ti) 
  newD$ESRD <- 0
  newD$YEARS_SURVIVAL <- tmax
  print(tmax)
  predSurv <-
    predict(
      JMmodel,
      newdata = newD,
      process = "event",
      times = seq(tmax, tmax + years_pred, length.out = 51),
      return_newdata = TRUE
    )
  predLong <-
    predict(
      JMmodel,
      newdata = newD,
      return_newdata = TRUE,
      times = seq(tmax, tmax + years_pred, length.out = 51) 
    ) 
  plot.predict_jm(predLong, predSurv,outcomes = c(2,1),
                  fun_long =  list(exp,identity),
                  ylab_long = c("eGFR", "PRO24H"),
                  #pos_ylab_long = c(0.1, 4),
                  cex_ylab_long = 1.6,
                  ylim_long_outcome_range=F,
                  tmax = tmax,years_pred= years_pred,
                  #fill_CI_event = "#FFEB77",
                  ylab_event = "Cumulative ESKD Risk"
                  
  )
  #return(predSurv)
}

#x = dynamicplotIgAN(JMmodel = jointFit_clinical5 ,newD = du,years_pred = 10)


###########################
###########################

plot.predict_jm <- function (x, x2 = NULL, subject = 1, outcomes = 1,
                             fun_long = NULL, fun_event = NULL,
                             CI_long = TRUE, CI_event = TRUE,
                             xlab = "Follow-up Time", ylab_long = NULL,
                             ylab_event = "Cumulative Risk", main = "",
                             lwd_long = 2, lwd_event = 2.5,
                             ylim_long_outcome_range = TRUE,
                             col_line_long = "#0000FF",
                             col_line_event = c("#FF0000", "#03BF3D", "#8000FF"),
                             pch_points = 16, col_points = "blue", cex_points = 1,
                             fill_CI_long = "#0000FF4D",
                             fill_CI_event = c("#FF00004D", "#03BF3D4D", "#8000FF4D"),
                             cex_xlab = 1, cex_ylab_long = 1, cex_ylab_event = 1.5,
                             cex_main = 1, cex_axis = 1, col_axis = "black",
                             pos_ylab_long = c(0.1, 2, 0.08), bg = "white",
                             tmax = tmax,years_pred = years_pred,
                             ...) {
  process_x <- attr(x, "process")
  pred_Long <- if (process_x == "longitudinal") x
  pred_Event <- if (process_x == "event") x
  if (!is.null(x2)) {
    process_x2 <- attr(x2, "process")
    if (process_x2 == "longitudinal" && is.null(pred_Long)) pred_Long <- x2
    if (process_x2 == "event" && is.null(pred_Event)) pred_Event <- x2
  }
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  resp_vars <- attr(x, "resp_vars")
  ranges <- attr(x, "ranges")
  last_times <- attr(x, "last_times")
  y <- attr(x, "y")
  times_y <- attr(x, "times_y")
  id <- attr(x, "id")
  if (!is.null(pred_Long)) {
    test1 <- is.data.frame(pred_Long)
    test2 <- is.list(pred_Long) && length(pred_Long) == 2L && is.data.frame(pred_Long[[1]])
    if (!test1 && !test2) {
      stop("you must use predict.jm(..., return_newdata = TRUE)")
    }
    if (test2) {
      pred_Long <- rbind(pred_Long[[1L]], pred_Long[[2L]])
    }
  }
  if (!is.null(pred_Event) && !is.data.frame(pred_Event)) {
    stop("you must use predict.jm(..., return_newdata = TRUE)")
  }
  unq_id <- if (!is.null(pred_Long)) unique(pred_Long[[id_var]])
  if (!is.null(pred_Event)) unq_id <- unique(c(pred_Event[[id_var]], unq_id))
  if (length(subject) > 1L) {
    stop("'subject' must be of length 1.")
  }
  if (!subject %in% unq_id && subject > length(unq_id)) {
    stop("not valid input for 'subject'.")
  }
  subj <- if (subject %in% unq_id) subject else unq_id[subject]
  subj_ind <- match(subj, unq_id)
  if (!is.null(pred_Long)) {
    pred_Long <- pred_Long[pred_Long[[id_var]] == subj, ]
    if (!is.null(pred_Event)) {
      pred_Long <- pred_Long[pred_Long[[time_var]] <= last_times[subj_ind], ]
    }
    if (!nrow(pred_Long)) {
      stop("no available measurements before the last time.")
    }
    pos_outcomes <- grep("pred_", names(pred_Long), fixed = TRUE)
    n_outcomes <- length(outcomes)
    if (n_outcomes > length(pos_outcomes)) {
      stop("the length of 'outcomes' is greater than the number of ",
           "outcomes in the dataset.")
    }
    if (any(outcomes > length(pos_outcomes))) {
      stop("not valid entries in 'outcome'.")
    }
    if (!is.null(pred_Event) && n_outcomes > 3) {
      warning("when 'pred_Event' is not null max three outcomes are allowed in the plot.")
      n_outcomes <- 3
      outcomes <- rep_len(outcomes, length.out = 3L)
    }
    if (is.null(fun_long)) {
      fun_long <- rep(list(function (x) x), n_outcomes)
    } else {
      if (is.function(fun_long)) fun_long <- rep(list(fun_long), n_outcomes)
      if (is.list(fun_long) && (length(fun_long) != n_outcomes ||
                                !all(sapply(fun_long, is.function)))) {
        stop("'fun_long' needs to be a function or a list of functions.")
      }
    }
    col_line_long <- rep(col_line_long, length.out = n_outcomes)
    pch_points <- rep(pch_points, length.out = n_outcomes)
    col_points <- rep(col_points, length.out = n_outcomes)
    cex_points <- rep(cex_points, length.out = n_outcomes)
    fill_CI_long <- rep(fill_CI_long, length.out = n_outcomes)
  }
  if (!is.null(pred_Event)) {
    pred_Event <- pred_Event[pred_Event[[id_var]] == subj, ]
    if (is.null(fun_event) || !is.function(fun_event)) {
      fun_event <- function (x) x
    }
  }
  if (is.null(ylab_long)) {
    ylab_long <- resp_vars
  }
  xlim <- NULL
  if (!is.null(pred_Long)) xlim <- range(xlim, pred_Long[[time_var]])
  if (!is.null(pred_Event)) xlim <- range(xlim, pred_Event[[time_var]])
  plot_long_i <- function (outcome, add_xlab = FALSE, box = TRUE,
                           cex_axis = cex_axis) {
    ind <- pos_outcomes[outcome]
    outcome_i <- match(outcome, outcomes)
    f <- fun_long[[outcome_i]]
    preds <- f(pred_Long[[ind]])
    low <- f(pred_Long[[ind + 1]])
    upp <- f(pred_Long[[ind + 2]])
    times <- pred_Long[[time_var]]
    na_preds <- is.na(preds)
    preds <- preds[!na_preds]
    low <- low[!na_preds]
    upp <- upp[!na_preds]
    times <- times[!na_preds]
    ry <- range(preds, low, upp, na.rm = TRUE)
    ry <- range(ry[1L] * 0.6, ry[2L] * 1.2) # <---
    rx <- range(times, na.rm = TRUE)
    y_lim <- if (ylim_long_outcome_range) {
      range(f(ranges[[outcome]]), ry, na.rm = TRUE)
    } else {
      ry
    }
    plot(rx, ry, type = "n", xaxt = "n", bty = if (box) "o" else "n",
         xlab = if (add_xlab) xlab  else "", xlim = xlim, col.axis = col_axis,
         ylim = y_lim, ylab = ylab_long[outcome],
         cex.lab = cex_ylab_long, cex.axis = cex_axis, col.lab = col_axis,
         col.axis = col_axis)
    if (!add_xlab) {
      axis(1, c(-5, last_times[subj_ind]), labels = c("", ""), tcl = 0,
           cex.axis = cex_axis, col = col_axis, col.axis = col_axis,
           col.ticks = col_axis)
    }
    if (CI_long) {
      polygon(c(times, rev(times)), c(low, rev(upp)), border = NA,
              col = fill_CI_long[outcome_i])
    }
    y_i <- f(c(y[[outcome]]))
    times_y_i <- times_y[[outcome]]
    id_i <- id[[outcome]]
    points(times_y_i[id_i == subj_ind], y_i[id_i == subj_ind],
           pch = pch_points[outcome_i], cex = cex_points[outcome_i],
           col = col_points[outcome_i])
    lines(times, preds, lwd = lwd_long, col = col_line_long[outcome_i])
    abline(v = last_times[subj_ind] + 0.01, lty = 3, col = col_axis)
  }
  plot_event <- function (box = FALSE, axis_side = 4, cex_axis = cex_axis) {
    ind <- grep("pred_", names(pred_Event), fixed = TRUE)
    preds <- fun_event(pred_Event[[ind]])
    low <- fun_event(pred_Event[[ind + 1]])
    upp <- fun_event(pred_Event[[ind + 2]])
    strata <- pred_Event[["_strata"]]
    if (is.null(strata)) strata <- rep(1, length(preds))
    unq_strata <- sort(unique(strata))
    col_line_event <- rep(col_line_event, length.out = length(unq_strata))
    fill_CI_event <- rep(fill_CI_event, length.out = length(unq_strata))
    times <- pred_Event[[time_var]] 
    ry <- sort(fun_event(c(0, 1)))
    rx <- range(times, na.rm = TRUE)
    plot(rx, ry, type = "n", xlab = "", ylab = "", xlim = xlim,
         axes = FALSE, col.axis = col_axis, col.lab = col_axis, ylim = ry)
    if (box) box(col = col_axis)
    axis(axis_side, cex.axis = cex_axis, col = col_axis,
         col.ticks = col_axis, col.axis = col_axis)
    for (i in seq_along(unq_strata)) {
      ind_str <- strata == unq_strata[i]
      if (CI_event) {
        polygon(c(times[ind_str], rev(times[ind_str])),
                c(low[ind_str], rev(upp[ind_str])), border = NA,
                col = fill_CI_event[i])
      }
      lines(times[ind_str], preds[ind_str], lwd = lwd_event,
            col = col_line_event[i])
    }
  }
  if (is.null(pred_Event)) {
    for (i in seq_along(outcomes)) {
      plot_long_i(outcomes[i], TRUE, cex_axis = cex_axis)
      title(main = main, cex = cex_main)
      axis(1, cex.axis = cex_axis, col = col_axis,
           col.ticks = col_axis, col.axis = col_axis)
    }
  }
  if (is.null(pred_Long)) {
    plot_event(box = TRUE, 2, cex_axis = cex_axis)
    title(xlab = xlab, cex = cex_xlab)
    title(ylab = ylab_event, cex = cex_ylab_event)
    title(main = main, cex = cex_main)
    abline(v = last_times[subj_ind] + 0.01, lty = 3)
    axis(1, cex.axis = cex_axis, col = col_axis, col.ticks = col_axis,
         col.axis = col_axis)
  }
  if (!is.null(pred_Long) && !is.null(pred_Event)) {
    if (n_outcomes == 1) {
      # n_outcomes == 1
      op <- par(mar = c(4,4,3,4), mgp = c(2, 0.4, 0), tcl = -0.3, bg = bg)
      plot_long_i(outcomes[1L], cex_axis = cex_axis)
      axis(1, cex.axis = cex_axis, col.ticks = col_axis, col = col_axis,
           col.axis = col_axis)
      title(xlab = xlab, cex = cex_xlab, col = col_axis)
      par(new = TRUE)
      plot_event(cex_axis = cex_axis)
      mtext(ylab_event, 4, 1.5, cex = cex_ylab_event, col = col_axis)
      par(op)
      mtext(main, 3, 1.5, cex = cex_main, col = col_axis)
    } else if (n_outcomes == 2) {
      # n_outcomes == 2
      op <- par(mfrow = c(2, 1), oma = c(4,4,3,4), mar = c(0, 0, 0, 0),
                mgp = c(2, 0.4, 0), tcl = -0.3, bg = bg)
      pp <- par("usr")[1] + pos_ylab_long * diff(par("usr")[1:2])
      plot_long_i(outcomes[1L], box = FALSE, cex_axis = cex_axis)
      mtext(ylab_long[outcomes[1L]], 2, 1.5, at = pp[1],
            cex = cex_ylab_long * 0.66, col = col_axis)
      plot_long_i(outcomes[2L], box = FALSE, cex_axis = cex_axis)
      mtext(ylab_long[outcomes[2L]], 2, 1.5, at = pp[2],
            cex = cex_ylab_long * 0.66, col = col_axis)
      axis(1, cex.axis = cex_axis, col.ticks = col_axis, col = col_axis,
           col.axis = col_axis)
      
      mtext(xlab, side = 1, line = 1.5, outer = TRUE,
            cex = cex_xlab, col = col_axis)
      par(op)
      op <- par(new = TRUE, oma = c(4,4,3,4), mar = c(0, 0, 0, 0),
                mgp = c(2, 0.4, 0), tcl = -0.3, cex = 0.9)
      plot_event(box = TRUE, cex_axis = 0.66 * cex_axis)
      mtext(ylab_event, 4, 1.5, cex = cex_ylab_event, col = col_axis)
      #abline(lty=2,v = tmax + 5)
      abline(lty=2,h = 0.5)
      axis(3, at= tmax+seq.int(0,years_pred), labels = seq.int(0,years_pred),
           cex.axis = 1.5*cex_axis, col.ticks = col_axis, col = col_axis, col.axis = col_axis)
      par(op)
      mtext(main, 3, 1.5, cex = cex_main, col = col_axis)
      
    } else {
      # n_outcomes == 3
      op <- par(mfrow = c(3, 1), oma = c(4,4,3,4), mar = c(0, 0, 0, 0),
                mgp = c(2, 0.4, 0), tcl = -0.3, bg = bg)
      pp <- par("usr")[1] + pos_ylab_long * diff(par("usr")[1:2])
      plot_long_i(outcomes[1L], box = FALSE, cex_axis = cex_axis)
      mtext(ylab_long[outcomes[1L]], 2, 1.5, at = pp[1],
            cex = cex_ylab_long * 0.66, col = col_axis)
      plot_long_i(outcomes[2L], box = FALSE, cex_axis = cex_axis)
      mtext(ylab_long[outcomes[2L]], 2, 1.5, at = pp[2],
            cex = cex_ylab_long * 0.66, col = col_axis)
      plot_long_i(outcomes[3L], box = FALSE, cex_axis = cex_axis)
      mtext(ylab_long[outcomes[3L]], 2, 1.5, at = pp[3],
            cex = cex_ylab_long * 0.66, col = col_axis)
      axis(1, cex.axis = cex_axis, col = col_axis, col.ticks = col_axis,
           col.axis = col_axis)
      mtext(xlab, side = 1, line = 1.5, outer = TRUE, cex = cex_xlab,
            col = col_axis)
      box("inner", col = col_axis)
      par(op)
      op <- par(new = TRUE, oma = 0.6525 * c(4,4,3,4), mar = c(0, 0, 0, 0),
                mgp = c(2, 0.4, 0), tcl = -0.3, cex = 0.66)
      plot_event(cex_axis = cex_axis)
      mtext(ylab_event, 4, 1.5, cex = cex_ylab_event, col = col_axis)
      par(op)
      mtext(main, 3, 1.5, cex = cex_main, col = col_axis)
    }
  }
  invisible()
}

# x = dynamicplotIgAN(JMmodel = jointFit_clinical5 ,newD = ND,years_pred = 8)
