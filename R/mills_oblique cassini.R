lat_ref <- 38
lon_ref <- -75
psi <- 45 * (pi / 180)
theta_ref <- lat_ref * (pi / 180)
phi_ref <- lon_ref * (pi / 180)


rho <- r <- 1
lat <- 40
lon <- -69

theta <- lat * (pi / 180)
phi <- lon * (pi / 180)


# From Mills 2014; doi: 10.1117/12.2060760

# SPHERICAL-CARTESIAN TRANSFORMS
## fwd proj in matrix form; rho = radius, theta = lat, phi = lon, Eq.2
r_vec <- function(rho, theta, phi) {
  matrix(
    c(
      rho * cos(theta) * cos(phi),
      rho * cos(theta) * sin(phi),
      sin(theta)
    ),
    dimnames = list(c("x", "y", "z"))
  )
}

rv <- r_vec(rho, theta, phi)
rv

## inv proj in matrix form; x/y/z coords on a sphere with radius rho, Eq. 3
# rho_vec = [ sqrt(x^2 + y^2 + z^2), asin(z/rho), atan(y/x) ]
rho_vec <- function(x, y, z, rho = 1) {
  matrix(
    c(
      sqrt(x^2 + y^2 + z^2),
      asin(z / rho),
      atan(y / x)
    ),
    dimnames = list(c('rho', 'theta', 'phi'))
  )
}

rho_vec(rv[1], rv[2], rv[3], rho)


## ROTATION MATRICES
R_x <- function(alpha) {
  matrix(
    c(
      1,
      0,
      0,
      0,
      cos(alpha),
      -sin(alpha),
      0,
      sin(alpha),
      cos(alpha)
    ),
    nrow = 3
  )
}

R_y <- function(alpha) {
  matrix(
    c(
      cos(alpha),
      0,
      sin(alpha),
      0,
      1,
      0,
      -sin(alpha),
      0,
      cos(alpha)
    ),
    nrow = 3
  )
}

R_z <- function(alpha) {
  matrix(
    c(
      cos(alpha),
      -sin(alpha),
      0,
      sin(alpha),
      cos(alpha),
      0,
      0,
      0,
      1
    ),
    nrow = 3
  )
}


# FORWARD OBLIQUE CASSINI
rho_triple_prime <- function(
  theta,
  phi,
  psi,
  theta_ref,
  phi_ref,
  rho = 6378137,
  s = 1
) {
  ## Step 1
  rho_prime <- matrix(c(rho, theta, phi - phi_ref))

  ## Step 2
  r_prime <- r_vec(rho_prime[1], rho_prime[2], rho_prime[3])
  # can also be written as r_prime <- do.call(r_vec, as.list(rho_prime))

  ## Step 3
  r_double_prime <- R_y(-theta_ref) %*% r_prime

  ## Step 4
  r_triple_prime <- R_x(-psi) %*% r_double_prime

  ## Step 5
  rho_triple_prime <- rho_vec(
    r_triple_prime[1],
    r_triple_prime[1],
    r_triple_prime[1],
    rho
  )

  rho_triple_prime
  ## Step 6
}


get_destination_point <- function(
  lat1,
  lon1,
  azimuth,
  distance_m,
  R = 6371000
) {
  # 1. Convert inputs to radians
  phi1 <- lat1 * pi / 180
  lam1 <- lon1 * pi / 180
  alpha <- azimuth * pi / 180
  d_r <- distance_m / R # Angular distance

  # 2. Spherical Law of Cosines for destination latitude
  phi2 <- asin(sin(phi1) * cos(d_r) + cos(phi1) * sin(d_r) * cos(alpha))

  # 3. Calculate destination longitude
  lam2 <- lam1 +
    atan2(sin(alpha) * sin(d_r) * cos(phi1), cos(d_r) - sin(phi1) * sin(phi2))

  # 4. Convert back to degrees and normalize longitude
  lat2_deg <- phi2 * 180 / pi
  lon2_deg <- ((lam2 * 180 / pi + 180) %% 360) - 180

  return(data.frame(lat = lat2_deg, lon = lon2_deg))
}

# Example: 100km away from 38N, 75W at 45 degrees azimuth
dest <- get_destination_point(38, -75, 45, 100000)
