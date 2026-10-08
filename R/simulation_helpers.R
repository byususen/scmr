# Data-generation helpers adapted from the supplied article script.
# All scenario parameters are explicit through control.

inside_irregular_domain <- function(coords, control) {
    control$domain_exclusion_s1_weight * coords[, 1]^2 + control$domain_exclusion_s2_weight * coords[,
        2]^2 > control$domain_exclusion_radius^2
}

sample_domain <- function(n, seed = NULL, control) {
    if (!is.null(seed))
        set.seed(seed)
    out <- matrix(NA_real_, nrow = n, ncol = 2)
    colnames(out) <- c("s1", "s2")
    got <- 0L
    while (got < n) {
        m <- max(1000L, 2L * (n - got))
        cand <- cbind(stats::runif(m, control$domain_s1_range[1], control$domain_s1_range[2]), stats::runif(m,
            control$domain_s2_range[1], control$domain_s2_range[2]))
        cand <- cand[inside_irregular_domain(cand, control = control), , drop = FALSE]
        take <- min(n - got, nrow(cand))
        if (take > 0) {
            out[(got + 1):(got + take), ] <- cand[seq_len(take), , drop = FALSE]
            got <- got + take
        }
    }
    out
}

region_bounds_6 <- function(control) {
    b <- control$region_s2_breaks
    data.frame(region = control$region_levels6, s1_min = c(rep(control$domain_s1_range[1], 3), rep(control$region_s1_split,
        3)), s1_max = c(rep(control$region_s1_split, 3), rep(control$domain_s1_range[2], 3)), s2_min = rep(b[1:3],
        2), s2_max = rep(b[2:4], 2), stringsAsFactors = FALSE)
}

region_label_6 <- function(coords, control) {
    s1 <- coords[, 1]
    s2 <- coords[, 2]
    j <- ifelse(s1 <= control$region_s1_split, 1L, 2L)
    k <- cut(s2, breaks = control$region_s2_breaks, include.lowest = TRUE, labels = FALSE)
    paste0("R", j, k)
}

region_index_6 <- function(coords, control) {
    s1 <- coords[, 1]
    s2 <- coords[, 2]
    j <- ifelse(s1 <= control$region_s1_split, 0L, 1L)
    k <- cut(s2, breaks = control$region_s2_breaks, include.lowest = TRUE, labels = FALSE) - 1L
    list(j = j, k = k)
}

allocate_counts_exact <- function(n, probs, control) {
    probs <- as.numeric(probs)
    probs <- probs/sum(probs)
    raw_counts <- n * probs
    counts <- floor(raw_counts)
    remainder <- n - sum(counts)
    if (remainder > 0) {
        frac <- raw_counts - counts
        add_idx <- order(frac, decreasing = TRUE)[seq_len(remainder)]
        counts[add_idx] <- counts[add_idx] + 1L
    }
    counts
}

sample_points_from_region <- function(n_g, bound_row, control) {
    if (n_g <= 0) {
        out <- matrix(numeric(0), nrow = 0, ncol = 2)
        colnames(out) <- c("s1", "s2")
        return(out)
    }
    out <- matrix(NA_real_, nrow = n_g, ncol = 2)
    colnames(out) <- c("s1", "s2")
    got <- 0L
    while (got < n_g) {
        m <- max(1000L, 3L * (n_g - got))
        cand <- cbind(stats::runif(m, bound_row$s1_min, bound_row$s1_max), stats::runif(m, bound_row$s2_min,
            bound_row$s2_max))
        cand <- cand[inside_irregular_domain(cand, control = control), , drop = FALSE]
        take <- min(n_g - got, nrow(cand))
        if (take > 0) {
            out[(got + 1):(got + take), ] <- cand[seq_len(take), , drop = FALSE]
            got <- got + take
        }
    }
    out
}

sample_domain_fixed_region_ratio <- function(n, region_probs, seed = NULL, control) {
    if (!is.null(seed))
        set.seed(seed)
    if (is.null(names(region_probs))) {
        names(region_probs) <- control$region_levels6
    }
    region_probs <- region_probs[control$region_levels6]
    region_probs <- region_probs/sum(region_probs)
    counts <- allocate_counts_exact(n, region_probs, control = control)
    names(counts) <- control$region_levels6
    bounds <- region_bounds_6(control = control)
    coords_list <- list()
    region_list <- list()
    for (rg in control$region_levels6) {
        b <- bounds[bounds$region == rg, , drop = FALSE]
        coords_g <- sample_points_from_region(n_g = counts[rg], bound_row = b, control = control)
        coords_list[[rg]] <- coords_g
        region_list[[rg]] <- rep(rg, nrow(coords_g))
    }
    coords <- do.call(rbind, coords_list)
    region <- unlist(region_list, use.names = FALSE)
    idx <- sample(seq_len(nrow(coords)))
    coords <- coords[idx, , drop = FALSE]
    region <- region[idx]
    list(coords = coords, region = region, target_probs = region_probs, target_counts = counts, realized_counts = table(factor(region,
        levels = control$region_levels6)))
}

choose_cluster_member_probs <- function(scenario_name, control) {
    if (scenario_name == "clustered_balanced") {
        return(control$cluster_member_probs_balanced)
    }
    if (scenario_name == "clustered_imbalanced") {
        return(control$cluster_member_probs_imbalanced)
    }
    return(NULL)
}

exp_cov_matrix <- function(coords, range_param, nugget = control$gp_nugget, kernel_power = 1, control) {
    D <- as.matrix(stats::dist(coords))
    Sigma <- exp(-(D/range_param)^kernel_power)
    diag(Sigma) <- diag(Sigma) + nugget
    Sigma
}

sim_gp <- function(coords, range_param, tau2 = 1, seed = NULL, kernel_power = 1, control) {
    if (!is.null(seed))
        set.seed(seed)
    Sigma <- tau2 * exp_cov_matrix(coords = coords, range_param = range_param, kernel_power = kernel_power,
        control = control)
    as.numeric(MASS::mvrnorm(n = 1, mu = rep(0, nrow(coords)), Sigma = Sigma))
}

simulate_predictors <- function(coords, p = 20, active = 1:5, rho_active = 0.5, rho_inactive = 0.1, eta = 0.6,
    predictor_eta_scale = 0.3, kernel_power = 2, seed = NULL, control) {
    if (!is.null(seed))
        set.seed(seed)
    n <- nrow(coords)
    X <- matrix(NA_real_, nrow = n, ncol = p)
    eta_eff <- eta * predictor_eta_scale
    if (!is.finite(eta_eff) || eta_eff <= 0) {
        stop("eta_eff must be positive.")
    }
    common_active <- sim_gp(coords = coords, range_param = eta_eff, tau2 = control$predictor_gp_tau2,
        kernel_power = kernel_power, control = control)
    common_inactive <- sim_gp(coords = coords, range_param = eta_eff, tau2 = control$predictor_gp_tau2,
        kernel_power = kernel_power, control = control)
    active_set <- sort(unique(active))
    inactive_set <- setdiff(seq_len(p), active_set)
    for (m in active_set) {
        zm <- sim_gp(coords = coords, range_param = eta_eff, tau2 = control$predictor_gp_tau2, kernel_power = kernel_power,
            control = control)
        X[, m] <- sqrt(rho_active) * common_active + sqrt(1 - rho_active) * zm
    }
    for (m in inactive_set) {
        zm <- sim_gp(coords = coords, range_param = eta_eff, tau2 = control$predictor_gp_tau2, kernel_power = kernel_power,
            control = control)
        X[, m] <- sqrt(rho_inactive) * common_inactive + sqrt(1 - rho_inactive) * zm
    }
    X <- as.matrix(scale(X))
    X[!is.finite(X)] <- 0
    colnames(X) <- paste0("x", seq_len(p))
    X
}

resize_active_pattern <- function(base_values, n_active, control) {
    n_active <- as.integer(n_active)
    if (n_active <= 0)
        return(numeric(0))
    base_values <- as.numeric(base_values)
    if (n_active <= length(base_values))
        return(base_values[seq_len(n_active)])
    out <- rep(base_values, length.out = n_active)
    scale_down <- 1/sqrt(seq_len(n_active)/length(base_values))
    scale_down[seq_len(length(base_values))] <- 1
    out * scale_down
}

simulate_beta_global <- function(coords, p = 20, active = 1:5, class_levels, control) {
    n <- nrow(coords)
    C <- length(class_levels)
    beta <- array(0, dim = c(n, C, p + 1), dimnames = list(NULL, class_levels, c("(Intercept)", paste0("x",
        seq_len(p)))))
    a_base <- resize_active_pattern(control$beta_global_class1_base, length(active), control = control)
    b_base <- resize_active_pattern(control$beta_global_class2_base, length(active), control = control)
    beta[, 1, 1] <- control$beta_global_intercept_class1
    beta[, 2, 1] <- control$beta_global_intercept_class2
    beta[, 3, 1] <- -(beta[, 1, 1] + beta[, 2, 1])
    for (jj in seq_along(active)) {
        m <- active[jj]
        beta[, 1, m + 1] <- a_base[jj]
        beta[, 2, m + 1] <- b_base[jj]
        beta[, 3, m + 1] <- -(beta[, 1, m + 1] + beta[, 2, m + 1])
    }
    beta
}

simulate_beta_clustered <- function(coords, p = 20, active = 1:5, class_levels, control) {
    n <- nrow(coords)
    C <- length(class_levels)
    beta <- array(0, dim = c(n, C, p + 1), dimnames = list(NULL, class_levels, c("(Intercept)", paste0("x",
        seq_len(p)))))
    idx <- region_index_6(coords, control = control)
    g1_vals <- control$beta_cluster_g1_values
    g2_vals <- control$beta_cluster_g2_values
    a_base <- resize_active_pattern(control$beta_cluster_class1_base, length(active), control = control)
    b_base <- resize_active_pattern(control$beta_cluster_class2_base, length(active), control = control)
    for (i in seq_len(n)) {
        g1 <- g1_vals[idx$j[i] + 1]
        g2 <- g2_vals[idx$k[i] + 1]
        q1 <- g1 + g2 + control$beta_cluster_q1_offset
        q2 <- g1 - g2 + control$beta_cluster_q2_offset
        beta[i, 1, 1] <- control$beta_cluster_intercept_scale_class1 * q1
        beta[i, 2, 1] <- control$beta_cluster_intercept_scale_class2 * q2
        beta[i, 3, 1] <- -(beta[i, 1, 1] + beta[i, 2, 1])
        for (jj in seq_along(active)) {
            m <- active[jj]
            beta[i, 1, m + 1] <- a_base[jj] * q1
            beta[i, 2, m + 1] <- b_base[jj] * q2
            beta[i, 3, m + 1] <- -(beta[i, 1, m + 1] + beta[i, 2, m + 1])
        }
    }
    beta
}

simulate_beta_smooth <- function(coords, p = 20, active = 1:5, seed = NULL, class_levels, psi_intercept = c(0.7,
    0.9), psi_active = c(0.5, 0.7, 0.9, 1.1, 1.3), beta_tau2_intercept = 0.8, beta_tau2_active = 0.8,
    center_gp = TRUE, control) {
    if (!is.null(seed))
        set.seed(seed)
    n <- nrow(coords)
    C <- length(class_levels)
    beta <- array(0, dim = c(n, C, p + 1), dimnames = list(NULL, class_levels, c("(Intercept)", paste0("x",
        seq_len(p)))))
    center_if_needed <- function(z) {
        if (isTRUE(center_gp))
            z - mean(z, na.rm = TRUE)
        else z
    }
    gp_int1 <- center_if_needed(sim_gp(coords = coords, range_param = psi_intercept[1], tau2 = beta_tau2_intercept,
        control = control))
    gp_int2 <- center_if_needed(sim_gp(coords = coords, range_param = psi_intercept[2], tau2 = beta_tau2_intercept,
        control = control))
    beta[, 1, 1] <- control$beta_global_intercept_class1 + gp_int1
    beta[, 2, 1] <- control$beta_global_intercept_class2 + gp_int2
    beta[, 3, 1] <- -(beta[, 1, 1] + beta[, 2, 1])
    a_base <- resize_active_pattern(control$beta_global_class1_base, length(active), control = control)
    b_base <- resize_active_pattern(control$beta_global_class2_base, length(active), control = control)
    for (jj in seq_along(active)) {
        m <- active[jj]
        gp1 <- center_if_needed(sim_gp(coords = coords, range_param = psi_active[jj], tau2 = beta_tau2_active,
            control = control))
        gp2 <- center_if_needed(sim_gp(coords = coords, range_param = psi_active[jj], tau2 = beta_tau2_active,
            control = control))
        beta[, 1, m + 1] <- a_base[jj] + gp1
        beta[, 2, m + 1] <- b_base[jj] + gp2
        beta[, 3, m + 1] <- -(beta[, 1, m + 1] + beta[, 2, m + 1])
    }
    beta
}

softmax_from_beta <- function(X, beta, control) {
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    C <- dim(beta)[2]
    eta <- matrix(0, nrow = n, ncol = C)
    for (cc in seq_len(C)) {
        beta_c <- matrix(beta[, cc, -1, drop = TRUE], nrow = n, ncol = p)
        eta[, cc] <- beta[, cc, 1] + rowSums(X * beta_c)
    }
    P <- softmax_rows(eta)
    colnames(P) <- dimnames(beta)[[2]]
    P
}

sample_multiclass <- function(P, seed = NULL, control) {
    if (!is.null(seed))
        set.seed(seed)
    classes <- colnames(P)
    apply(P, 1, function(pr) {
        sample(classes, size = 1, prob = pr)
    })
}

get_class_balance_shift <- function(class_balance, K, control) {
    if (class_balance == "balanced") {
        return(rep(0, K))
    }
    if (class_balance == "imbalanced") {
        if (K == 3) {
            return(control$class_imbalance_shift_three)
        }
        else {
            z <- seq(control$class_imbalance_shift_range_general[1], control$class_imbalance_shift_range_general[2],
                length.out = K)
            return(z - mean(z))
        }
    }
    stop("Unknown class balance scenario: ", class_balance)
}

make_simulation_dataset <- function(n_obs, n_new, p, active, eta, scenario, class_balance, seed, class_levels,
    smooth_setup = NULL, predictor_eta_scale = 0.3, predictor_kernel_power = 2, predictor_rho_active = 0.5,
    predictor_rho_inactive = 0.1, cluster_member_probs = NULL, heterogeneity_strength = 1, control) {
    set.seed(seed)
    if (!scenario %in% c("global", "clustered_balanced", "clustered_imbalanced", "smooth", "clustered_irregular")) {
        stop("Unknown scenario: ", scenario)
    }
    if (!is.finite(heterogeneity_strength) || heterogeneity_strength < 0) {
        stop("heterogeneity_strength must be finite and nonnegative.")
    }
    if (is.null(cluster_member_probs)) {
        coords_all <- sample_domain(n_obs + n_new, seed = seed + 1, control = control)
        region_all <- region_label_6(coords_all, control = control)
        obs_idx <- seq_len(n_obs)
        new_idx <- if (n_new > 0)
            n_obs + seq_len(n_new)
        else integer(0)
        observed_region_counts <- table(factor(region_all[obs_idx], levels = control$region_levels6))
        new_region_counts <- table(factor(region_all[new_idx], levels = control$region_levels6))
        target_probs <- rep(NA_real_, length(control$region_levels6))
        names(target_probs) <- control$region_levels6
    }
    else {
        obs_coord_obj <- sample_domain_fixed_region_ratio(n = n_obs, region_probs = cluster_member_probs,
            seed = seed + 1, control = control)
        new_coord_obj <- sample_domain_fixed_region_ratio(n = n_new, region_probs = cluster_member_probs,
            seed = seed + 11, control = control)
        coords_all <- rbind(obs_coord_obj$coords, new_coord_obj$coords)
        region_all <- c(obs_coord_obj$region, new_coord_obj$region)
        obs_idx <- seq_len(n_obs)
        new_idx <- if (n_new > 0)
            n_obs + seq_len(n_new)
        else integer(0)
        observed_region_counts <- obs_coord_obj$realized_counts
        new_region_counts <- new_coord_obj$realized_counts
        target_probs <- cluster_member_probs[control$region_levels6]
        target_probs <- target_probs/sum(target_probs)
    }
    true_cluster_all <- switch(scenario, global = rep("T1", n_obs + n_new), clustered_balanced = as.character(region_all),
        clustered_imbalanced = as.character(region_all), smooth = rep(NA_character_, n_obs + n_new),
        clustered_irregular = paste0("Q", irregular_regime_index(coords_all, control) + 1L))
    X_all <- simulate_predictors(coords = coords_all, p = p, active = active, rho_active = predictor_rho_active,
        rho_inactive = predictor_rho_inactive, eta = eta, predictor_eta_scale = predictor_eta_scale,
        kernel_power = predictor_kernel_power, seed = seed + 2, control = control)
    if (is.null(smooth_setup))
        smooth_setup <- list()
    smooth_setup_use <- smooth_setup
    smooth_setup_use$control <- control
    if (scenario == "smooth") {
        n_act <- length(active)
        smooth_setup_use$psi_active <- if (n_act <= 1) {
            rep(control$smooth_psi_active_min, max(n_act, 1))
        }
        else {
            seq(control$smooth_psi_active_min, control$smooth_psi_active_max, length.out = n_act)
        }
    }
    beta_all <- switch(scenario, global = simulate_beta_global(coords = coords_all, p = p, active = active,
        class_levels = class_levels, control = control), clustered_balanced = {
        beta_global_ref <- simulate_beta_global(coords = coords_all, p = p, active = active, class_levels = class_levels,
            control = control)
        beta_cluster_raw <- simulate_beta_clustered(coords = coords_all, p = p, active = active, class_levels = class_levels,
            control = control)
        beta_global_ref + heterogeneity_strength * (beta_cluster_raw - beta_global_ref)
    }, clustered_imbalanced = {
        beta_global_ref <- simulate_beta_global(coords = coords_all, p = p, active = active, class_levels = class_levels,
            control = control)
        beta_cluster_raw <- simulate_beta_clustered(coords = coords_all, p = p, active = active, class_levels = class_levels,
            control = control)
        beta_global_ref + heterogeneity_strength * (beta_cluster_raw - beta_global_ref)
    }, clustered_irregular = {
        beta_global_ref <- simulate_beta_global(coords = coords_all, p = p, active = active, class_levels = class_levels,
            control = control)
        beta_regime_raw <- simulate_beta_irregular(coords = coords_all, p = p, active = active, class_levels = class_levels,
            control = control)
        beta_global_ref + heterogeneity_strength * (beta_regime_raw - beta_global_ref)
    }, smooth = do.call(simulate_beta_smooth, c(list(coords = coords_all, p = p, active = active, seed = seed +
        3, class_levels = class_levels), smooth_setup_use)))
    class_shift <- get_class_balance_shift(class_balance = class_balance, K = length(class_levels), control = control)
    for (cc in seq_along(class_levels)) {
        beta_all[, cc, 1] <- beta_all[, cc, 1] + class_shift[cc]
    }
    P_all <- softmax_from_beta(X_all, beta_all, control = control)
    Y_all <- factor(sample_multiclass(P_all, seed = seed + 4, control = control), levels = class_levels)
    list(observed = list(X = X_all[obs_idx, , drop = FALSE], y = Y_all[obs_idx], coords = coords_all[obs_idx,
        , drop = FALSE], beta = beta_all[obs_idx, , , drop = FALSE], prob = P_all[obs_idx, , drop = FALSE],
        region = region_all[obs_idx], true_cluster = true_cluster_all[obs_idx], unit = paste0("U", obs_idx)),
        new = list(X = X_all[new_idx, , drop = FALSE], y = Y_all[new_idx], coords = coords_all[new_idx,
            , drop = FALSE], beta = beta_all[new_idx, , , drop = FALSE], prob = P_all[new_idx, , drop = FALSE],
            region = region_all[new_idx], true_cluster = true_cluster_all[new_idx], unit = paste0("N",
                seq_along(new_idx))), class_shift = class_shift, cluster_member_probs = cluster_member_probs,
        target_region_probs = target_probs, observed_region_counts = observed_region_counts, new_region_counts = new_region_counts,
        heterogeneity_strength = heterogeneity_strength, n_active = length(active), n_inactive = p -
            length(active))
}

# Latin-square regimes on a grid of tiles: regime = (row + column) mod 3.
irregular_regime_index <- function(coords, control) {
    nt <- control$irregular_tiles
    col <- floor((coords[, 1] - control$domain_s1_range[1]) / diff(control$domain_s1_range) * nt[1])
    row <- floor((coords[, 2] - control$domain_s2_range[1]) / diff(control$domain_s2_range) * nt[2])
    col <- pmin(pmax(col, 0), nt[1] - 1)
    row <- pmin(pmax(row, 0), nt[2] - 1)
    as.integer((row + col) %% 3)
}

simulate_beta_irregular <- function(coords, p = 20, active = 1:5, class_levels, control) {
    n <- nrow(coords)
    C <- length(class_levels)
    beta <- array(0, dim = c(n, C, p + 1), dimnames = list(NULL, class_levels, c("(Intercept)", paste0("x",
        seq_len(p)))))
    regime <- irregular_regime_index(coords, control)
    q <- matrix(control$irregular_regime_q, ncol = 2, byrow = TRUE)
    a_base <- resize_active_pattern(control$beta_cluster_class1_base, length(active), control = control)
    b_base <- resize_active_pattern(control$beta_cluster_class2_base, length(active), control = control)
    for (i in seq_len(n)) {
        q1 <- q[regime[i] + 1L, 1]
        q2 <- q[regime[i] + 1L, 2]
        beta[i, 1, 1] <- control$beta_cluster_intercept_scale_class1 * q1
        beta[i, 2, 1] <- control$beta_cluster_intercept_scale_class2 * q2
        beta[i, 3, 1] <- -(beta[i, 1, 1] + beta[i, 2, 1])
        for (jj in seq_along(active)) {
            m <- active[jj]
            beta[i, 1, m + 1] <- a_base[jj] * q1
            beta[i, 2, m + 1] <- b_base[jj] * q2
            beta[i, 3, m + 1] <- -(beta[i, 1, m + 1] + beta[i, 2, m + 1])
        }
    }
    beta
}
