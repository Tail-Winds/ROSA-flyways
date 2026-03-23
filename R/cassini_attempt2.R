oblique_cassini <- function(lat, lon, lat_ref, lon_ref, psi, radius = 6378137) {
  # lat, lon:        Input vectors (degrees)
  # lat_ref, lon_ref: Center of projection (degrees)
  # psi:             Azimuth of great circle (degrees)
  # radius:          Earth radius in meters (default WGS84 semi-major)

  d2r <- pi / 180

  # 1. Convert to Radians
  theta <- lat * d2r
  phi <- lon * d2r
  t_ref <- lat_ref * d2r
  p_ref <- lon_ref * d2r
  psi_r <- psi * d2r

  # 2. Cartesian Conversion (on Unit Sphere, rho=1)
  # We work on the unit sphere for rotations, then scale at the end.
  # r_vec = [x, y, z] rows
  # Note: Standard mapping x=cos(lat)cos(lon), y=cos(lat)sin(lon), z=sin(lat)
  r_vec <- rbind(
    cos(theta) * cos(phi),
    cos(theta) * sin(phi),
    sin(theta)
  )

  # 3. Define Rotation Matrices

  # Rot_z(-phi_ref)
  Rz <- matrix(
    c(cos(-p_ref), -sin(-p_ref), 0, sin(-p_ref), cos(-p_ref), 0, 0, 0, 1),
    nrow = 3,
    byrow = TRUE
  )

  # Rot_y(theta_ref)
  Ry <- matrix(
    c(cos(t_ref), 0, sin(t_ref), 0, 1, 0, -sin(t_ref), 0, cos(t_ref)),
    nrow = 3,
    byrow = TRUE
  )

  # Rot_x(-psi)
  Rx <- matrix(
    c(1, 0, 0, 0, cos(-psi_r), -sin(-psi_r), 0, sin(-psi_r), cos(-psi_r)),
    nrow = 3,
    byrow = TRUE
  )

  # 4. Apply Rotations
  # Sequence: Z (long shift) -> Y (lat shift) -> X (track align)
  r_final <- Rx %*% Ry %*% Rz %*% r_vec

  # 5. Extract Rotated Coordinates
  x_prime <- r_final[1, ]
  y_prime <- r_final[2, ]
  z_prime <- r_final[3, ]

  # 6. Convert to Spherical Coordinates (Radians)
  # theta_rot: Latitude in rotated system (distance from great circle)
  # phi_rot:   Longitude in rotated system (distance along great circle)

  # Clamp z to [-1, 1] to avoid NaNs from floating point noise
  theta_rot <- asin(pmax(pmin(z_prime, 1), -1))
  phi_rot <- atan2(y_prime, x_prime)

  # 7. Apply Earth Radius Scaling
  # Distance = Angle (radians) * Radius

  x_map <- radius * phi_rot
  y_map <- radius * theta_rot

  return(data.frame(x = x_map, y = y_map))
}

oblique_cassini(45, -78.1, 45, -78, 45)

inverse_oblique_cassini <- function(
  x,
  y,
  lat_ref,
  lon_ref,
  psi,
  radius = 6378137
) {
  # x, y:              Input map coordinates (meters)
  # lat_ref, lon_ref:  Center of projection (degrees)
  # psi:               Azimuth of great circle (degrees)
  # radius:            Earth radius (meters)

  d2r <- pi / 180
  r2d <- 180 / pi

  # 1. Convert Map Coordinates to Radians (Rotated Spherical)
  # x = radius * phi_rot  -> phi_rot = x / radius
  # y = radius * theta_rot -> theta_rot = y / radius

  phi_rot <- x / radius
  theta_rot <- y / radius

  # 2. Convert Rotated Spherical to Cartesian Vector (r''')
  # r''' = [x''', y''', z''']
  # Standard mapping: x = cos(lat)cos(lon), y = cos(lat)sin(lon), z = sin(lat)

  x_tprime <- cos(theta_rot) * cos(phi_rot)
  y_tprime <- cos(theta_rot) * sin(phi_rot)
  z_tprime <- sin(theta_rot)

  r_vec_tprime <- rbind(x_tprime, y_tprime, z_tprime)

  # 3. Define Inverse Rotation Matrices
  # The inverse of a rotation matrix is its transpose.
  # Or simply rotate by the negative angle.

  # Psi, Theta_ref, Phi_ref in radians
  psi_rad <- psi * d2r
  theta_ref <- lat_ref * d2r
  phi_ref <- lon_ref * d2r

  # Inverse Rot_x(-psi) -> Rotate X by +psi
  Rx_inv <- matrix(
    c(1, 0, 0, 0, cos(psi_rad), -sin(psi_rad), 0, sin(psi_rad), cos(psi_rad)),
    nrow = 3,
    byrow = TRUE
  )

  # Inverse Rot_y(theta_ref) -> Rotate Y by -theta_ref
  Ry_inv <- matrix(
    c(
      cos(-theta_ref),
      0,
      sin(-theta_ref),
      0,
      1,
      0,
      -sin(-theta_ref),
      0,
      cos(-theta_ref)
    ),
    nrow = 3,
    byrow = TRUE
  )

  # Inverse Rot_z(-phi_ref) -> Rotate Z by +phi_ref
  Rz_inv <- matrix(
    c(cos(phi_ref), -sin(phi_ref), 0, sin(phi_ref), cos(phi_ref), 0, 0, 0, 1),
    nrow = 3,
    byrow = TRUE
  )

  # 4. Apply Inverse Transformation Chain
  # Sequence is reversed: Rx_inv -> Ry_inv -> Rz_inv
  # r_orig = Rz_inv . Ry_inv . Rx_inv . r'''

  r_orig <- Rz_inv %*% Ry_inv %*% Rx_inv %*% r_vec_tprime

  # 5. Cartesian to Geographic (Lat/Lon)
  x_final <- r_orig[1, ]
  y_final <- r_orig[2, ]
  z_final <- r_orig[3, ]

  # Latitude = asin(z)
  lat_rad <- asin(pmax(pmin(z_final, 1), -1))

  # Longitude = atan2(y, x)
  lon_rad <- atan2(y_final, x_final)

  # Convert to Degrees
  return(data.frame(
    lat = lat_rad * r2d,
    lon = lon_rad * r2d
  ))
}


inverse_oblique_cassini(-5562.54, 5569.407, 45, -78, 45)
