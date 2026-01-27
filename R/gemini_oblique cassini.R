# Oblique Cassini Projection for a Sphere
# phi, lambda: Lat/Lon of point to project (radians)
# phi_c, lambda_c: Lat/Lon of the central origin (radians)
# alpha: Azimuth of the central line at the origin (radians, clockwise from North)
# R: Radius of the sphere (e.g., 6371000 meters)

project_oblique_cassini <- function(
  phi,
  lambda,
  phi_c,
  lambda_c,
  alpha,
  R = 6371000
) {
  # 1. Calculate the Pole (phi_p, lambda_p) of the oblique great circle
  phi_p <- asin(cos(phi_c) * sin(alpha))
  lambda_p <- atan2(cos(alpha), -sin(phi_c) * sin(alpha)) + lambda_c

  # 2. Transform to Oblique Coordinates (phi_prime, lambda_prime)
  # These are the unsimplified spherical rotation equations
  sin_phi_prime <- sin(phi) *
    sin(phi_p) +
    cos(phi) * cos(phi_p) * sin(lambda - lambda_p)
  phi_prime <- asin(sin_phi_prime)

  # Components for lambda_prime
  y_comp <- cos(phi) * cos(lambda - lambda_p)
  x_comp <- cos(phi_p) *
    sin(phi) -
    sin(phi_p) * cos(phi) * sin(lambda - lambda_p)
  lambda_prime <- atan2(y_comp, x_comp)

  # 3. Apply Cassini Forward Formulas
  # B is the sine of the angular distance from the central line
  B <- cos(phi_prime) * sin(lambda_prime)

  x <- R * asin(B)
  y <- R * atan2(tan(phi_prime), cos(lambda_prime))

  # 4. Shift y so that the origin (phi_c, lambda_c) is at y = 0
  # We calculate the y-value of the center point in this projection space
  # (Alternatively, you can just return x and y as raw coordinates)

  return(data.frame(x = x, y = y))
}

# --- Example Usage ---
deg2rad <- function(deg) deg * pi / 180

# Define a path starting in London (51.5N, 0W) heading North-East (45 degrees)
center_lat <- deg2rad(51.5)
center_lon <- deg2rad(0)
azimuth <- deg2rad(45)

# Point to project (e.g., Paris: 48.85N, 2.35E)
test_lat <- deg2rad(48.85)
test_lon <- deg2rad(2.35)

coords <- project_oblique_cassini(
  test_lat,
  test_lon,
  center_lat,
  center_lon,
  azimuth
)
print(coords)


# Spherical Cassini Projection
# phi, lambda: Coordinates of the point (in radians)
# phi_0, lambda_0: Coordinates of the origin (in radians)
# R: Radius of the sphere (default is WGS84 mean radius in meters)

project_cassini <- function(phi, lambda, phi_0, lambda_0, R = 6371007) {
  # 1. Calculate the intermediate value B
  # B is the sine of the angular distance from the central meridian
  B <- cos(phi) * sin(lambda - lambda_0)

  # 2. Forward Formulas (Snyder 13-1 and 13-2)
  # x is the perpendicular distance to the central meridian
  x <- R * asin(B)

  # y is the distance along the central meridian from the origin latitude
  y <- R * (atan2(tan(phi), cos(lambda - lambda_0)) - phi_0)

  return(data.frame(x = x, y = y))
}

# --- Example Usage ---
deg2rad <- function(deg) deg * pi / 180

# Origin: Greenwich (0, 0)
lat_0 <- deg2rad(0)
lon_0 <- deg2rad(0)

# Point: Paris (48.8566 N, 2.3522 E)
p_lat <- deg2rad(48.8566)
p_lon <- deg2rad(2.3522)

coords <- project_cassini(p_lat, p_lon, lat_0, lon_0)
print(coords)
# Output x is approx 172 km (east of central meridian)
# Output y is approx 5433 km (north of equator)

#
project_azimuth_cassini <- function(
  phi,
  lambda,
  phi_c,
  lambda_c,
  alpha,
  R = 6371000
) {
  # 1. Pole Calculation
  phi_p <- asin(cos(phi_c) * sin(alpha)) #Eq 9-7
  lambda_p <- atan2(cos(alpha), -sin(phi_c) * sin(alpha)) + lambda_c # Eq 9-8

  # Helper to get raw coordinates
  get_raw <- function(p, l) {
    # 2. Spherical Rotation
    s_phi_prime <- sin(p) * sin(phi_p) - cos(p) * cos(phi_p) * sin(l - lambda_p)
    phi_prime <- asin(s_phi_prime)

    y_c <- cos(p) * cos(l - lambda_p)
    x_c <- sin(p) * cos(phi_p) + cos(p) * sin(phi_p) * sin(l - lambda_p)
    lambda_prime <- atan2(y_c, x_c)

    # 3. Cassini Formulas
    raw_x <- R * asin(cos(phi_prime) * sin(lambda_prime))
    raw_y <- R * atan2(tan(phi_prime), cos(lambda_prime))
    return(c(raw_x, raw_y))
  }

  # Center the origin
  offset <- get_raw(phi_c, lambda_c)
  results <- mapply(get_raw, phi, lambda)

  return(data.frame(x = results[1, ] - offset[1], y = results[2, ] - offset[2]))
}


# Inverse Oblique Cassini Projection (Spherical)
inverse_oblique_cassini <- function(x, y, phi_c, lambda_c, alpha, R = 6371000) {
  # 1. Re-calculate the Pole (must match the forward step)
  phi_p <- asin(cos(phi_c) * sin(alpha))
  lambda_p <- atan2(cos(alpha), -sin(phi_c) * sin(alpha)) + lambda_c

  # 2. Re-calculate the y_offset used in the forward step
  # We need the raw y value of the center point to "un-shift" our input y
  get_raw_y_offset <- function(p, l) {
    s_phi_prime <- sin(p) * sin(phi_p) - cos(p) * cos(phi_p) * sin(l - lambda_p)
    phi_prime <- asin(s_phi_prime)
    y_c <- cos(p) * cos(l - lambda_p)
    x_c <- sin(p) * cos(phi_p) + cos(p) * sin(phi_p) * sin(l - lambda_p)
    lambda_prime <- atan2(y_c, x_c)
    return(R * atan2(tan(phi_prime), cos(lambda_prime)))
  }

  y_offset <- get_raw_y_offset(phi_c, lambda_c)

  # 3. Process the points
  inverse_calc <- function(xi, yi) {
    y_raw <- yi + y_offset

    # Oblique components
    D <- y_raw / R
    x_rad <- xi / R

    phi_prime <- asin(sin(D) * cos(x_rad))
    lambda_prime <- atan2(tan(x_rad), cos(D))

    # 4. Inverse Spherical Rotation
    phi <- asin(
      sin(phi_prime) *
        sin(phi_p) +
        cos(phi_prime) * cos(phi_p) * cos(lambda_prime)
    )

    y_term <- cos(phi_prime) * sin(lambda_prime)
    x_term <- sin(phi_prime) *
      cos(phi_p) -
      cos(phi_prime) * sin(phi_p) * cos(lambda_prime)
    lambda <- lambda_p + atan2(y_term, x_term)

    return(c(lat = phi * 180 / pi, lon = lambda * 180 / pi))
  }

  results <- mapply(inverse_calc, x, y)
  return(as.data.frame(t(results)))
}

phi <- 39
lambda <- -74.5
phi_c <- 38 * pi / 180
lambda_c <- -75 * pi / 180
alpha <- 45 * pi / 180
proj <- project_azimuth_cassini(
  39 * pi / 180,
  -74.5 * pi / 180,
  phi_c,
  lambda_c,
  alpha
)
inverse_oblique_cassini(proj[1], proj[2], phi_c, lambda_c, alpha)


####
create_oblique_cassini <- function(lat_c, lon_c, azimuth, R = 6371000) {
  # 1. Convert to Radians
  phi_c <- lat_c * pi / 180
  lam_c <- lon_c * pi / 180
  alpha <- azimuth * pi / 180

  # 2. Define the Rotation Pole (Snyder's standard oblique method)
  # This pole defines the new "Equator" which is your path
  phi_p <- asin(cos(phi_c) * sin(alpha))
  lam_p <- atan2(cos(alpha), -sin(phi_c) * sin(alpha)) + lam_c

  # 3. Calculate the Longitude Shift (to force origin x to 0)
  # We find the oblique longitude of the center point
  l_obs_c <- atan2(
    cos(phi_c) * sin(lam_c - lam_p),
    sin(phi_c) * cos(phi_p) - cos(phi_c) * sin(phi_p) * cos(lam_c - lam_p)
  )

  # --- Forward Projection ---
  project <- function(lat, lon) {
    p <- lat * pi / 180
    l <- lon * pi / 180

    # Geographic -> Oblique Spherical
    sin_p_obs <- sin(p) * sin(phi_p) + cos(p) * cos(phi_p) * cos(l - lam_p)
    p_obs <- asin(sin_p_obs)

    # Oblique longitude relative to the pole
    l_obs_raw <- atan2(
      cos(p) * sin(l - lam_p),
      sin(p) * cos(phi_p) - cos(p) * sin(phi_p) * cos(l - lam_p)
    )

    # Shift longitude so center is 0
    l_obs <- l_obs_raw - l_obs_c

    # Cassini Plane Equations (Equidistant)
    x <- R * asin(cos(p_obs) * sin(l_obs))
    y <- R * atan2(tan(p_obs), cos(l_obs))

    return(data.frame(x = x, y = y))
  }

  # --- Inverse Projection ---
  inverse <- function(x, y) {
    # Cassini -> Oblique Spherical
    p_obs <- asin(sin(y / R) * cos(x / R))
    l_obs_shifted <- atan2(tan(x / R), cos(y / R))

    # Un-shift longitude
    l_obs <- l_obs_shifted + l_obs_c

    # Oblique -> Geographic (Full Inverse Rotation)
    # Using the Law of Cosines/Sines for spherical triangles
    phi <- asin(sin(p_obs) * sin(phi_p) + cos(p_obs) * cos(phi_p) * cos(l_obs))

    d_lam <- atan2(
      cos(p_obs) * sin(l_obs),
      sin(p_obs) * cos(phi_p) - cos(p_obs) * sin(phi_p) * cos(l_obs)
    )
    lam <- d_lam + lam_p

    # Normalize
    lon_deg <- ((lam * 180 / pi + 180) %% 360) - 180
    return(data.frame(lat = phi * 180 / pi, lon = lon_deg))
  }

  return(list(project = project, inverse = inverse))
}

my_map <- create_oblique_cassini(38, -75, 45)

# 1. Test Origin
origin <- my_map$project(38, -75)
print(origin)

# 2. Test Round Trip
fwd <- my_map$project(40, -74)
back <- my_map$inverse(fwd$x, fwd$y)
print(back)


#########
create_oblique_cassini <- function(lat_c, lon_c, azimuth, R = 6371000) {
  # 1. Convert to Radians
  phi_c <- lat_c * pi / 180
  lam_c <- lon_c * pi / 180
  alpha <- azimuth * pi / 180

  # 2. Define the Rotation Pole
  phi_p <- asin(cos(phi_c) * sin(alpha))
  lam_p <- atan2(cos(alpha), -sin(phi_c) * sin(alpha)) + lam_c

  # 3. Calculate Center Longitude Offset
  l_obs_c <- atan2(
    cos(phi_c) * sin(lam_c - lam_p),
    sin(phi_c) * cos(phi_p) - cos(phi_c) * sin(phi_p) * cos(lam_c - lam_p)
  )

  # --- Forward Projection ---
  project <- function(lat, lon) {
    p <- lat * pi / 180
    l <- lon * pi / 180

    sin_p_obs <- sin(p) * sin(phi_p) + cos(p) * cos(phi_p) * cos(l - lam_p)
    p_obs <- asin(sin_p_obs)

    l_obs_raw <- atan2(
      cos(p) * sin(l - lam_p),
      sin(p) * cos(phi_p) - cos(p) * sin(phi_p) * cos(l - lam_p)
    )
    l_obs <- l_obs_raw - l_obs_c

    x <- R * asin(cos(p_obs) * sin(l_obs))
    y <- R * atan2(tan(p_obs), cos(l_obs))

    # Clean up near-zero values
    x[abs(x) < 1e-9] <- 0
    y[abs(y) < 1e-9] <- 0

    return(data.frame(x = x, y = y))
  }

  # --- Inverse Projection ---
  inverse <- function(x, y) {
    p_obs <- asin(sin(y / R) * cos(x / R))
    l_obs <- atan2(tan(x / R), cos(y / R)) + l_obs_c

    phi <- asin(sin(p_obs) * sin(phi_p) + cos(p_obs) * cos(phi_p) * cos(l_obs))
    d_lam <- atan2(
      cos(p_obs) * sin(l_obs),
      sin(p_obs) * cos(phi_p) - cos(p_obs) * sin(phi_p) * cos(l_obs)
    )
    lam <- d_lam + lam_p

    return(data.frame(
      lat = phi * 180 / pi,
      lon = ((lam * 180 / pi + 180) %% 360) - 180
    ))
  }

  # --- Grid Plotting Helper ---
  plot_grid <- function(lat_range, lon_range, step = 1) {
    lats <- seq(lat_range[1], lat_range[2], by = step)
    lons <- seq(lon_range[1], lon_range[2], by = step)

    grid_lines <- list()

    # Latitude lines (parallels)
    for (la in lats) {
      seq_lon <- seq(lon_range[1], lon_range[2], length.out = 50)
      df <- project(rep(la, 50), seq_lon)
      df$id <- paste0("lat_", la)
      grid_lines[[length(grid_lines) + 1]] <- df
    }

    # Longitude lines (meridians)
    for (lo in lons) {
      seq_lat <- seq(lat_range[1], lat_range[2], length.out = 50)
      df <- project(seq_lat, rep(lo, 50))
      df$id <- paste0("lon_", lo)
      grid_lines[[length(grid_lines) + 1]] <- df
    }

    all_lines <- do.call(rbind, grid_lines)
    library(ggplot2)
    ggplot(all_lines, aes(x / 1000, y / 1000, group = id)) +
      geom_path(color = "grey80") +
      geom_point(aes(x = 0, y = 0), color = "red", size = 3) +
      coord_fixed() +
      labs(title = "Oblique Cassini Grid", x = "X (km)", y = "Y (km)") +
      theme_minimal()
  }

  return(list(project = project, inverse = inverse, plot_grid = plot_grid))
}
my_map <- create_oblique_cassini(38, -75, 45)

my_map$project(38, -75)
my_map$project(60, -75)
my_map$inverse(1707980, 1772902)
# Plot a 10x10 degree grid around the origin
my_map$plot_grid(c(33, 43), c(-80, -70), step = 1)


library(ggplot2)

# 1. Initialize your specific map
# Center: 38N, 75W | Path Azimuth: 45 degrees
my_map <- create_oblique_cassini(lat_c = 38, lon_c = -75, azimuth = 45)

# 2. Function to create a geographic circle (Indicatrix)
generate_indicatrix <- function(lat, lon, radius_km = 20) {
  # Create 60 points around a circle
  theta <- seq(0, 2 * pi, length.out = 60)

  # Approximate conversion: 1 degree lat approx 111.3km
  # Adjusting lon radius by cos(lat) to keep it a circle on the sphere
  d_lat <- radius_km / 111.3
  d_lon <- radius_km / (111.3 * cos(lat * pi / 180))

  circle_lat <- lat + d_lat * cos(theta)
  circle_lon <- lon + d_lon * sin(theta)

  # Project these points into our Cassini x,y
  proj <- my_map$project(circle_lat, circle_lon)
  proj$id <- paste0(lat, "_", lon)
  return(proj)
}

# 3. Define a grid of points to place the circles
# We'll look at a 4x4 degree area around the center
lats <- seq(36, 45, by = 1)
lons <- seq(-77, -69, by = 1)
grid <- expand.grid(lat = lats, lon = lons)

# 4. Generate all ellipses
tissot_data <- do.call(
  rbind,
  lapply(1:nrow(grid), function(i) {
    generate_indicatrix(grid$lat[i], grid$lon[i], radius_km = 15)
  })
)

# 5. Plot the result
ggplot(tissot_data, aes(x = x / 1000, y = y / 1000, group = id)) +
  geom_polygon(fill = "steelblue", color = "blue", alpha = 0.4) +
  # Add the central path (x=0)
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  coord_fixed() +
  labs(
    title = "Tissot's Indicatrices: Oblique Cassini Projection",
    subtitle = "Distortion increases as vertical stretching away from the central red path",
    x = "Cross-track Distance (km)",
    y = "Along-track Distance (km)"
  ) +
  theme_minimal()


#####
# Create map centered at 38, -75 with azimuth 45
my_map <- create_oblique_cassini(38, -75, 52)

# We know the 'oblique' frame has the path along its equator.
# We will pick a point at y=200km, then move it exactly 50km perpendicular.

# Point A: On the path
pt_a_geo <- my_map$inverse(0, 200000)
# Point B: Exactly 50km perpendicular (x=50000)
pt_b_geo <- my_map$inverse(50000, 200000)

# Now, calculate the Great Circle distance on the sphere between these two Geo points
# Using the Haversine formula
dist_haversine <- function(lat1, lon1, lat2, lon2, R = 6371000) {
  phi1 <- lat1 * pi / 180
  phi2 <- lat2 * pi / 180
  dphi <- (lat2 - lat1) * pi / 180
  dlam <- (lon2 - lon1) * pi / 180
  a <- sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlam / 2)^2
  return(R * 2 * atan2(sqrt(a), sqrt(1 - a)))
}

true_dist <- dist_haversine(
  pt_a_geo$lat,
  pt_a_geo$lon,
  pt_b_geo$lat,
  pt_b_geo$lon
)

cat("Map X-Distance:", 50000, "meters\n")
cat("True Sphere Distance:", true_dist, "meters\n")
cat("Scale Factor (h):", 50000 / true_dist)

# check scale on main axis
# Point A: On the path
pt_a_geo <- my_map$inverse(0, 0)
# Point B: Exactly 50km perpon meridian (y=50000)
pt_b_geo <- my_map$inverse(0, 50000)
dist_haversine(pt_a_geo$lat, pt_a_geo$lon, pt_b_geo$lat, pt_b_geo$lon)
