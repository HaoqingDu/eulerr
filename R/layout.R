#' Skyline packing algorithm
#'
#' @param m a matrix of vertices for the rectangles
#'
#' @return A matrix with updated vertices for the rectangles.
#'
#' @keywords internal
skyline_pack <- function(m) {
  # TODO: Add rotation to boxes as well.
  # TODO: Port to c++

  n <- NCOL(m)
  w <- m[2L, ] - m[1L, ]
  h <- m[4L, ] - m[3L, ]
  sizes <- h * w

  # Add some padding for the rectangles
  #padding <- 1.6*sqrt(sum(sizes))*0.015
  padding <- sum(w, na.rm = TRUE) * 0.8 * 0.015

  # Pick a maximum bin width. Make sure the largest rectangle fits.
  #bin_w <- max(1.6*sqrt(sum(sizes)), w + padding)
  bin_w <- max(
    sum(w, na.rm = TRUE) * 0.8,
    sum(w[order(w, decreasing = TRUE)][1:2] + 2 * padding),
    w + padding
  )

  w <- w + padding
  h <- h + padding

  tol <- sqrt(.Machine$double.eps)

  m[] <- 0

  # Initialize the skyline
  skyline <- matrix(c(0, 0, bin_w, 0), ncol = 2)

  for (i in 1L:n) {
    # Order the points by y coordinate
    ord <- order(skyline[2L, ])

    # Start by examining the lowest rooftop on the skyline
    k <- 0L

    looking <- TRUE
    while (looking) {
      j <- which(ord == (1L + 2L * k) | ord == (2L + 2L * k))[1L]

      p1 <- ord[ord == j]
      p2 <- ord[ord == j + 1L]

      left <- skyline[2L, 1L:p1] - skyline[2L, p1] > tol
      right <- skyline[2L, p2:NCOL(skyline)] - skyline[2L, p2] > tol

      if (any(left)) {
        # There is a taller rooftop on the skyline to the left
        next_left <- which(left)[sum(left)]
      } else {
        next_left <- 1L
      }

      if (any(right)) {
        # There is a taller rooftop on the skyline to the right
        next_right <- which(right)[1L]
      } else {
        next_right <- NCOL(skyline)
      }

      if (w[i] <= skyline[1L, next_right] - skyline[1L, next_left]) {
        # Build a new building on the skyline
        m[1L, i] <- skyline[1L, next_left]
        m[2L, i] <- skyline[1L, next_left] + w[i]
        m[3L, i] <- skyline[2L, p1]
        m[4L, i] <- skyline[2L, p2] + h[i]

        l <- if (next_left == 1) 0 else 1

        skyline[2L, next_left + l] <- skyline[2L, p1] + h[i]

        newcols <-
          rbind(
            c(skyline[1L, next_left] + w[i], skyline[1L, next_left] + w[i]),
            c(skyline[2L, p2] + h[i], skyline[2L, p2])
          )

        skyline <- cbind(
          skyline[, seq(1L, next_left + l), drop = FALSE],
          newcols,
          skyline[, seq(p2, NCOL(skyline)), drop = FALSE]
        )

        # Check if there are any rooftops on the skyline beneath the new one
        underneath <- skyline[1L, ] > m[1L, i] & skyline[1L, ] < m[2L, i]

        if (any(underneath)) {
          # Drop down to the lowest level
          skyline[2L, which(underneath)[1L] - 1L] <-
            skyline[2L, which(underneath)[sum(underneath)]]
          skyline <- skyline[, !underneath, drop = FALSE]
        }

        looking <- FALSE
      } else {
        # Examine the next rooftop
        k <- k + 1L
      }
    }
  }
  m
}

#' Compress an Euler layout
#'
#' @param fpar an Euler layout fit with [euler()]
#' @param id the binary index of sets
#'
#' @return A modified fpar object.
#' @keywords internal
# compress_layout <- function(fpar, id, fit) {
#   # TODO: Port to c++
#   n <- NCOL(id)
# 
#   clusters <- matrix(NA, nrow = n, ncol = n)
# 
#   for (i in 1:n) {
#     for (j in 1:n) {
#       clusters[i, j] <- any(id[, i] & id[, j] & fit > 0)
#     }
#   }
# 
#   for (i in 1:n) {
#     for (j in 1:n) {
#       if (any(clusters[i, ] & clusters[j, ])) {
#         clusters[i, ] <- clusters[j, ] <- clusters[i, ] | clusters[j, ]
#       }
#     }
#   }
# 
#   unique_clusters <- unique(lapply(split(clusters, row(clusters)), which))
# 
#   # Drop clusters that contain no elements (usually shapes without area)
#   unique_clusters <- unique_clusters[lengths(unique_clusters) > 0L]
#   n_clusters <- length(unique_clusters)
# 
#   if (n_clusters > 0) {
#     bounds <- matrix(NA, ncol = n_clusters, nrow = 4L)
# 
#     for (i in seq_along(unique_clusters)) {
#       ii <- unique_clusters[[i]]
# 
#       h <- fpar[ii, 1]
#       k <- fpar[ii, 2]
#       a <- fpar[ii, 3]
#       b <- fpar[ii, 4]
#       phi <- fpar[ii, 5]
# 
#       # normalize rotation by setting rotation angle between two first
#       # ellipses to 0
# 
#       if (length(h) > 1) {
#         theta <- atan2(k[2] - k[1], h[2] - h[1])
# 
#         h0 <- cos(-theta)*(h - h[1]) - sin(-theta)*(k - k[1]) + h[1]
#         k0 <- sin(-theta)*(h - h[1]) + cos(-theta)*(k - k[1]) + k[1]
#         phi0 <- phi - theta
# 
#         h <- h0
#         k <- k0
#         phi <- phi0
# 
#         xc <- mean(range(h))
#         yc <- mean(range(k))
# 
#         # mirror across y axis if first shape is not at bottom
#         if ((k[1] > yc)) {
#           k <- yc - k
#           phi <- pi - phi
#         }
# 
#         # mirror across x axis if first set is not furthest to the left
#         if (h[1] > xc) {
#           h <- xc - h
#           phi <- pi - phi
#         }
#       }
# 
#       limits <- get_bounding_box(h, k, a, b, phi)
# 
#       fpar[ii, 1] <- h
#       fpar[ii, 2] <- k
#       fpar[ii, 3] <- a
#       fpar[ii, 4] <- b
#       fpar[ii, 5] <- phi
# 
#       bounds[1:2, i] <- limits$xlim
#       bounds[3:4, i] <- limits$ylim
#     }
# 
#     if (n_clusters > 1) {
#       # pack the bounding rectangles
#       # TODO: Fix occasional errors in computing the bounding boxes.
#       if (all(is.finite(bounds))) {
#         new_bounds <- skyline_pack(bounds)
#         for (i in seq_along(unique_clusters)) {
#           ii <- unique_clusters[[i]]
#           fpar[ii, 1] <- fpar[ii, 1] - (bounds[1, i] - new_bounds[1, i])
#           fpar[ii, 2] <- fpar[ii, 2] - (bounds[3, i] - new_bounds[3, i])
#         }
#       }
#     }
#   }
#   fpar
# }

compress_layout <- function(fpar, id, fit) {
  n <- NCOL(id)
  clusters <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      clusters[i, j] <- any(id[, i] & id[, j] & fit > 0)
    }
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (any(clusters[i, ] & clusters[j, ])) {
        clusters[i, ] <- clusters[j, ] <- clusters[i, ] | clusters[j, ]
      }
    }
  }
  
  unique_clusters <- unique(lapply(split(clusters, row(clusters)), which))
  unique_clusters <- unique_clusters[lengths(unique_clusters) > 0L]
  n_clusters <- length(unique_clusters)
  
  if (n_clusters > 0) {
    bounds <- matrix(NA, ncol = n_clusters, nrow = 4L)
    
    for (i in seq_along(unique_clusters)) {
      ii <- unique_clusters[[i]]
      
      h <- fpar[ii, 1]
      k <- fpar[ii, 2]
      a <- fpar[ii, 3]
      b <- fpar[ii, 4]
      phi <- fpar[ii, 5]
      
      if (length(h) > 1) {
        theta <- atan2(k[2] - k[1], h[2] - h[1])       

        h0 <- cos(-theta) * (h - h[1]) - sin(-theta) * (k - k[1]) + h[1]
        k0 <- sin(-theta) * (h - h[1]) + cos(-theta) * (k - k[1]) + k[1]

        phi0 <- phi - theta
        
        h <- h0
        k <- k0
        phi <- phi0
        
        xc <- mean(range(h))
        yc <- mean(range(k))
        
        if ((k[1] > yc)) {
          k <- yc - k
          phi <- pi - phi
        }
        
        if (h[1] > xc) {
          h <- xc - h
          phi <- pi - phi
        }
      }
      
      limits <- get_bounding_box(h, k, a, b, phi)
      
      fpar[ii, 1] <- h
      fpar[ii, 2] <- k
      fpar[ii, 3] <- a
      fpar[ii, 4] <- b
      fpar[ii, 5] <- phi
      
      bounds[1:2, i] <- limits$xlim
      bounds[3:4, i] <- limits$ylim
    }
    
    if (n_clusters > 1) {
      if (all(is.finite(bounds))) {
        new_bounds <- skyline_pack(bounds)
        for (i in seq_along(unique_clusters)) {
          ii <- unique_clusters[[i]]
          fpar[ii, 1] <- fpar[ii, 1] - (bounds[1, i] - new_bounds[1, i])
          fpar[ii, 2] <- fpar[ii, 2] - (bounds[3, i] - new_bounds[3, i])
        }
      }
    }
  }
  
  # === ADDITION: Rotate layout so h[1] == h[2] ===
  if (nrow(fpar) >= 2) {
    h1 <- fpar[1, 1]
    k1 <- fpar[1, 2]
    h2 <- fpar[2, 1]
    k2 <- fpar[2, 2]
    
    theta <- atan2(k2 - k1, h2 - h1)
    cos_theta <- cos(-theta)
    sin_theta <- sin(-theta)
    
    h_rot <- cos_theta * (fpar[, 1] - h1) - sin_theta * (fpar[, 2] - k1) + h1
    k_rot <- sin_theta * (fpar[, 1] - h1) + cos_theta * (fpar[, 2] - k1) + k1
    phi_rot <- fpar[, 5] - theta
    
    fpar[, 1] <- h_rot
    fpar[, 2] <- k_rot
    fpar[, 5] <- phi_rot
  }
  
  # === Mirror k[3] and/or k[4] if needed ===
  if (nrow(fpar) >= 4) {
    k1 <- fpar[1, 2]
    for (i in 3:4) {
      if (!is.na(fpar[i, 2])) {
        k_val <- fpar[i, 2]
        if ((i == 3 && k_val < k1) || (i == 4 && k_val >= k1)) {
          fpar[i, 2] <- 2 * k1 - k_val
          fpar[i, 5] <- pi - fpar[i, 5]
        }
      }
    }
  }
  
  # === Check angle at point 1 and spin point 3 if obtuse ===
  h1 <- fpar[1, 1]; k1 <- fpar[1, 2]
  h2 <- fpar[2, 1]; k2 <- fpar[2, 2]
  h3 <- fpar[3, 1]; k3 <- fpar[3, 2]
  
  d12_sq <- (h2 - h1)^2 + (k2 - k1)^2
  d13_sq <- (h3 - h1)^2 + (k3 - k1)^2
  d23_sq <- (h3 - h2)^2 + (k3 - k2)^2
  
  a2 <- fpar[2, 3]; b2 <- fpar[2, 4]
  a3 <- fpar[3, 3]; b3 <- fpar[3, 4]
  r2 <- mean(c(a2, b2))
  r3 <- mean(c(a3, b3))
  radii_sum_sq_32<- (r2 + r3)^2
  
  if (d23_sq > radii_sum_sq_32) {
    
    max_d <- max(d12_sq, d13_sq, d23_sq)
    sum_others <- d12_sq + d13_sq + d23_sq - max_d
    
    if (max_d > sum_others) {
      # Compute angle at 1 using dot product
      v1x <- h2 - h1; v1y <- k2 - k1
      v2x <- h3 - h1; v2y <- k3 - k1
      dot <- v1x * v2x + v1y * v2y
      norm1 <- sqrt(d12_sq)
      norm2 <- sqrt(d13_sq)
      cos_angle <- dot / (norm1 * norm2)
      angle1 <- acos(cos_angle)
      delta_angle <- angle1 - pi * 0.4  # Adjust to 90°
      
      cos_d <- cos(-delta_angle)
      sin_d <- sin(-delta_angle)
      h3_new <- cos_d * (h3 - h1) - sin_d * (k3 - k1) + h1
      k3_new <- sin_d * (h3 - h1) + cos_d * (k3 - k1) + k1
      phi3_new <- fpar[3, 5] - delta_angle
      
      # fpar[3, 1] <- h3_new
      # fpar[3, 2] <- k3_new
      # fpar[3, 5] <- phi3_new
      
      d23_new_sq <- (h3_new - h2)^2 + (k3_new - k2)^2
      
      if (d23_new_sq > radii_sum_sq_32) {
        fpar[3, 1] <- h3_new
        fpar[3, 2] <- k3_new
        fpar[3, 5] <- phi3_new
      } else {
        # Compute vector from 2 to original 3
        dx <- h3 - h2
        dy <- k3 - k2
        norm <- sqrt(dx^2 + dy^2)
        
        if (norm == 0) {
          dx <- 1  # Arbitrary direction to avoid zero division
          dy <- 0
          norm <- 1
        }
        
        dx <- dx / norm
        dy <- dy / norm
        
        # Set distance = r2 + r3
        h3_alt <- h2 + dx * (r2 + r3)
        k3_alt <- k2 + dy * (r2 + r3)
        
        # Preserve original angle for phi
        fpar[3, 1] <- h3_alt
        fpar[3, 2] <- k3_alt
        # Optionally update phi to match direction
        fpar[3, 5] <- atan2(dy, dx)
      }
      
    }
  }
  
  # === Adjust point 4 to top of point 2 ===
  h1 <- fpar[1, 1]; k1 <- fpar[1, 2]
  h2 <- fpar[2, 1]; k2 <- fpar[2, 2]
  h4 <- fpar[4, 1]; k4 <- fpar[4, 2]
  
  d12_sq <- (h2 - h1)^2 + (k2 - k1)^2
  d14_sq <- (h4 - h1)^2 + (k4 - k1)^2
  d24_sq <- (h4 - h2)^2 + (k4 - k2)^2
  
  a1 <- fpar[1, 3]; b1 <- fpar[1, 4]
  a4 <- fpar[4, 3]; b4 <- fpar[4, 4]
  r1 <- mean(c(a1, b1))
  r4 <- mean(c(a4, b4))
  radii_sum_sq_41 <- (r1 + r4)^2
  
  if (d14_sq > radii_sum_sq_41) {
    
    max_d <- max(d12_sq, d14_sq, d24_sq)
    sum_others <- d12_sq + d14_sq + d24_sq - max_d
    
    if (max_d > sum_others) {
      v1x <- h1 - h2; v1y <- k1 - k2
      v2x <- h4 - h2; v2y <- k4 - k2
      dot <- v1x * v2x + v1y * v2y
      norm1 <- sqrt(d12_sq)
      norm2 <- sqrt(d24_sq)
      cos_angle <- dot / (norm1 * norm2)
      angle2 <- acos(cos_angle)
      delta_angle <- angle2 - pi * 0.4
      
      cos_d <- cos(-delta_angle)
      sin_d <- sin(-delta_angle)
      h4_new <- cos_d * (h4 - h2) - sin_d * (k4 - k2) + h2
      k4_new <- sin_d * (h4 - h2) + cos_d * (k4 - k2) + k2
      phi4_new <- fpar[4, 5] - delta_angle
      
      d14_new_sq <- (h4_new - h1)^2 + (k4_new - k1)^2
      
      if (d14_new_sq > radii_sum_sq_41) {
        fpar[4, 1] <- h4_new
        fpar[4, 2] <- k4_new
        fpar[4, 5] <- phi4_new
      } else {
        # Compute vector from 1 to original 4
        dx <- h4 - h1
        dy <- k4 - k1
        norm <- sqrt(dx^2 + dy^2)
        
        if (norm == 0) {
          dx <- 0
          dy <- 1  # Arbitrary vertical direction
          norm <- 1
        }
        
        dx <- dx / norm
        dy <- dy / norm
        
        # Set distance = r1 + r4
        h4_alt <- h1 + dx * (r1 + r4)
        k4_alt <- k1 + dy * (r1 + r4)
        
        # Set position and optionally adjust phi to match direction
        fpar[4, 1] <- h4_alt
        fpar[4, 2] <- k4_alt
        fpar[4, 5] <- atan2(dy, dx)
      }
      
    }
  }
  
  fpar
}

#' Center ellipses
#'
#' @param pars a matrix or data.frame of x coordinates, y coordinates, minor
#'   radius (a) and major radius (b)
#'
#' @return A centered version of `pars`.
#' @keywords internal
center_layout <- function(pars) {
  void <- pars[, 3L] == 0 | pars[, 4L] == 0
  h <- pars[!void, 1L]
  k <- pars[!void, 2L]
  a <- pars[!void, 3L]
  b <- pars[!void, 4L]
  phi <- pars[!void, 5L]

  cphi <- cos(phi)
  sphi <- sin(phi)
  xlim <- range(c(h + a * cphi, h + b * cphi, h - a * cphi, h - b * cphi))
  ylim <- range(c(k + a * sphi, k + b * sphi, k - a * sphi, k - b * sphi))

  pars[!void, 1L] <- h + abs(xlim[1L] - xlim[2L]) / 2 - xlim[2L]
  pars[!void, 2L] <- k + abs(ylim[1L] - ylim[2L]) / 2 - ylim[2L]
  if (any(!void)) {
    pars[void, 1L] <- sum(pars[!void, 1L]) / sum(!void)
    pars[void, 2L] <- sum(pars[!void, 2L]) / sum(!void)
  }
  pars
}
