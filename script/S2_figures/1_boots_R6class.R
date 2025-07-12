# Bootstrapping Object 
# Newest version: (2025/04/12)

library(pbapply) # needed for progress bar in *apply function
library(R6) # needed for OOP
library(tidyverse)
library(broom)

BaseBootstrap <- R6Class("BaseBootstrap",
  public = list(
    # Data and model parameters
    data = NULL,
    y_var = NULL,
    x_vars = NULL,
    
    # Results storage
    boot_results = NULL,
    aggregated_data = NULL,
    x_range = NULL,
    
    # Initialize with data and model specifications
    initialize = function(data, y_var, x_vars){
      self$data <- data
      self$y_var <- y_var
      self$x_vars <- x_vars
    }
  ),
  
  private = list(
    create_sim_mx = function(x_max, x_min, tmax, tmin, order) {
      # order=c(p, q) the numbers of pairs and poly; ncol = 2q+p+1
      p <- order[1]; q <- order[2]
      r <- seq(x_min, x_max, 0.5)
      s <- (r - tmin)/(tmax - tmin)
      mx1 <- matrix(sapply(1:p, function(i) s^i), ncol=p, byrow=FALSE)
      mx2 <- matrix(sapply(1:q, function(i) c(sin(2*i*pi*s), cos(2*i*pi*s))), ncol=2*q, byrow=FALSE)
      ones <- matrix(1, nrow=length(s), ncol=1)
      cbind(ones, mx1, mx2)
    },
    
    fill_table = function(boot_results, sim_mx) {
      # Use stored x_range
      rounds <- length(boot_results)
      x_min <- self$x_range[1]
      x_max <- self$x_range[2]
      
      # Build an empty dataframe
      df <- data.frame(
        matrix(NA, nrow=nrow(sim_mx), ncol=(rounds+1))
      )
      colnames(df) <- c("temperature", paste0("round_", seq_len(rounds)))
      # Fill in temp
      df$temperature <- seq(x_min, x_max, 0.5)
      # Loop: every round has a set of fitted values
      for (rd in 1:rounds){
        coef_vec <- unname(boot_results[[rd]]) # strip names
        ghat_x <- sim_mx %*% coef_vec
        j <- rd + 1
        df[, j] <- ghat_x[, 1]
      }
      # Calculate quantiles
      df$CI_lb <- apply(df[, -1], 1, quantile, probs=0.025, na.rm=TRUE)
      df$fv_median <- apply(df[, -1], 1, median, na.rm=TRUE)
      df$CI_ub <- apply(df[, -1], 1, quantile, probs=0.975, na.rm=TRUE)
      
      return(df)
    }
  )
)

# TRF Class (inherits Base) -------------------------------------
TRF_Model <- R6Class("TRF_Model",
  inherit = BaseBootstrap,
  
  public = list(
    sim_matrix = NULL,
    order = NULL,
    temp_vars = NULL,
    # Create simulation matrix
    create_sim_matrix = function(x_range, t_range, order) {
      self$x_range <- x_range
      self$order <- order
      x_min <- x_range[1]; x_max <- x_range[2]
      t_min <- t_range[1]; t_max <- t_range[2]
      self$sim_matrix <- private$create_sim_mx(x_max, x_min, t_max, t_min, order)
      # Automatically detect temperature variables
      self$temp_vars <- grep("ts|sin|cos", self$x_vars, value = TRUE)
    },
    
    # Implement TRF Bootstrap
    run_bootstrap = function(rounds) {
      y <- self$data[[self$y_var]]
      X <- model.matrix(~ ., data = self$data[self$x_vars])
      rgr <- lm.fit(X, y) # faster than lm()
      u <- residuals(rgr)
      
      # Bootstrap loop - TRF
      self$boot_results <- pblapply(1:rounds, function(rd) {
        u_ <- sample(u, replace = TRUE)
        y_ <- X %*% rgr$coefficient + u_
        lm.fit(X, y_)$coefficients # faster
      })
    },
    
    aggregate_results = function() {
      # Convert to plain numeric vectors
      boot_results_subset <- lapply(self$boot_results, function(coefs) {
        unname(unlist(c(coefs["(Intercept)"], coefs[self$temp_vars])))
      })
      # Verify dimensions
      expected_length <- ncol(self$sim_matrix)
      if (!all(sapply(boot_results_subset, length) == expected_length)) {
        stop(paste("Coefficient length mismatch. Expected", expected_length, 
                   "but got", unique(sapply(boot_results_subset, length))))
      }
      self$aggregated_data <- private$fill_table(boot_results_subset, self$sim_matrix)
    },
    
    # TRF plotting: Add to existing layer
    plot = function(base_plot, save_path = NULL) {
      p <- base_plot+
        geom_line(data=self$aggregated_data, aes(x=temperature, y=fv_median))+
        geom_ribbon(data=self$aggregated_data, 
                    aes(x=temperature, ymin = CI_lb, ymax = CI_ub), 
                    fill="grey", alpha=0.5)
      if (!is.null(save_path)) ggsave(save_path, p)
      return(p)
    }
  )
)

# CTRF Class (inherits Base) -------------------------------------
CTRF_Model <- R6Class("CTRF_Model",
  inherit = BaseBootstrap,
  
  public = list(
    var_groups = NULL,  # List of variable subsets, list(vars1, vars2, ...)
    sim_matrices = NULL,  # store multiple simulation matrices
    orders = NULL,  # List of orders, list(c(2,2), c(1,1), ...)
    
    initialize = function(data, y_var, x_vars, var_groups, orders) {
      super$initialize(data, y_var, x_vars)
      self$var_groups <- var_groups
      self$orders <- orders
    },
    
    create_sim_matrices = function(x_range, t_range) {
      self$x_range <- x_range
      self$sim_matrices <- lapply(seq_along(self$var_groups), function(i) {
        order <- self$orders[[i]]
        private$create_sim_mx(x_range[2], x_range[1], t_range[2], t_range[1], order)
      })
    },
    
    # Implement CTRF Bootstrap
    run_bootstrap = function(rounds) {
      y <- self$data[[self$y_var]]
      X <- model.matrix(~ ., data = self$data[self$x_vars])
      rgr <- lm.fit(X, y)
      u <- residuals(rgr)
      
      # Bootstrap loop - CTRF
      self$boot_results <- pblapply(1:rounds, function(rd) {
        u_ <- sample(u, replace = TRUE)
        y_ <- X %*% rgr$coefficient + u_
        coefs <- lm.fit(X, y_)$coefficients
        
        # Return a list of coefficients by var_groups
        lapply(seq_along(self$var_groups), function(i) {
          vars <- self$var_groups[[i]]
          if (i == 1) { # for base TRF
            c(coefs["(Intercept)"], coefs[vars])
          } else { # for covariate groups
            coefs[vars]
          }
        })
      })
    },
    
    aggregate_results = function() {
      # Process each variable group separately
      self$aggregated_data <- lapply(seq_along(self$var_groups), function(i) {
        group_results <- lapply(self$boot_results, `[[`, i) 
        ## extraction operator [[ to extract the i-th element from each list in boot_results
        private$fill_table(group_results, self$sim_matrices[[i]]) # use group-specific matrix
      })
    },
    
    # CTRF plotting
    plot = function(save_dir, x_breaks=seq(-5, 35, 5), hrsuffix=NULL) {
      for (i in 2:length(self$aggregated_data)) {
        filename <- paste0("CTRF_", i, hrsuffix, ".png")
        p <- ggplot(self$aggregated_data[[i]], aes(temperature, fv_median)) +
          geom_line() +
          geom_hline(yintercept=0, color="plum", linetype="solid")+
          geom_ribbon(aes(ymin=CI_lb, ymax=CI_ub), fill="grey", alpha=0.5)+
          labs(x="Temperature (C)", y="Partial Effect")+
          scale_x_continuous(breaks = x_breaks)
        ggsave(file.path(save_dir, filename), p, width=5, height=0.618*5)
      }
    }
  )
)

#################################################################
# If automatic plotting of CTRF is unsatisfactory, do it manually. 
# Not for Base TRF
drawsave_CTRF <- function(df, w_name, v_breaks, suffix, save_dir, y_lim=NULL){
  filename <- paste0("CTRF_", suffix, ".png")
  plot_CTRF <- ggplot()+
    geom_line(data=df, aes(x=temperature, y=fv_median))+
    geom_ribbon(data=df, aes(x=temperature, ymin=CI_lb, ymax=CI_ub), fill="grey", alpha=0.5)+
    geom_hline(yintercept=0, color="plum", linetype="solid")+
    labs(x="Temperature (C)", y=paste0(w_name, " partial effect"))+
    scale_x_continuous(breaks = v_breaks)
  if (!is.null(y_lim)) {
    plot_CTRF + scale_y_continuous(limits=y_lim)
  }
  ggsave(file = file.path(save_dir, filename), width=5, height=0.618*5)
  invisible()
}

# Draw hourly plot
# This could be used as 'base_plot' in TRF_Model environment
draw_tempplot <- function(data, hour, y_var, y_lab, color_hour, x_lim, y_lim){
  df <- dplyr::filter(data, Hour==hour) # Drawing scatterplot doesn't need separate dt_hr beforehand
  plot <- ggplot()+ 
    geom_point(mapping = aes(x=df$temperature, y=df[[y_var]]), color=color_hour, alpha=0.3)+
    labs(x="Temperature (C)", y=y_lab)+
    scale_x_continuous(limits=x_lim)+ # Fixed limits make it easier to compare hours 
    scale_y_continuous(limits=y_lim)+
    theme_classic()
  return(plot)
}
