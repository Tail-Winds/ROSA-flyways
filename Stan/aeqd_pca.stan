functions {
  // Returns N x 2 matrix of projected coordinates
  matrix project_aeqd(vector lat, vector lon, real lat0, real lon0) {
    int N = rows(lat);
    matrix[N, 2] coords;
    for (n in 1:N) {
      real cos_c = sin(lat0) * sin(lat[n]) + cos(lat0) * cos(lat[n]) * cos(lon[n] - lon0);
      real c = acos(fmin(1.0, fmax(-1.0, cos_c)));
      real k_prime = (c < 1e-9) ? 1.0 : c / sin(c);
      
      coords[n, 1] = k_prime * cos(lat[n]) * sin(lon[n] - lon0);
      coords[n, 2] = k_prime * (cos(lat0) * sin(lat[n]) - sin(lat0) * cos(lat[n]) * cos(lon[n] - lon0));
    }
    return coords;
  }
}

data {
  int<lower=1> N;
  vector[N] lat_rad;
  vector[N] lon_rad;
}

parameters {
  unit_vector[3] center_vec; 
}

transformed parameters {
  real lat0 = asin(center_vec[3]);
  real lon0 = atan2(center_vec[2], center_vec[1]);
}

model {
  // Finding the Geodesic Median (optimal center)
  for (n in 1:N) {
    real cos_c = center_vec[1] * cos(lat_rad[n]) * cos(lon_rad[n]) + 
                 center_vec[2] * cos(lat_rad[n]) * sin(lon_rad[n]) + 
                 center_vec[3] * sin(lat_rad[n]);
    target += -acos(fmin(1.0, fmax(-1.0, cos_c)));
  }
}

generated quantities {
  matrix[N, 2] projected_points = project_aeqd(lat_rad, lon_rad, lat0, lon0);
  
  // 1. Center the projected coordinates for PCA
  vector[2] mu_proj;
  matrix[N, 2] centered_points;
  for (j in 1:2) {
    mu_proj[j] = mean(projected_points[, j]);
    centered_points[, j] = projected_points[, j] - mu_proj[j];
  }
  
  // 2. Compute the 2x2 Covariance Matrix
  matrix[2, 2] cov_mat = (centered_points' * centered_points) / (N - 1);
  
  // 3. Eigen-decomposition (PCA)
  // eigenvectors[, 1] is the 1st Principal Component
  complex_matrix[2, 2] eigvects_complex = eigenvectors(cov_mat);
  complex_vector[2] eigvals_complex = eigenvalues(cov_mat);

  // 4. Extract the real parts for use in standard real types
  matrix[2, 2] eigvects = get_real(eigvects_complex);
  vector[2] eigvals = get_real(eigvals_complex);

  // 5. Calculate the alignment angle (in degrees)
  // This tells you how much to rotate (+gamma) to level the data
  real pca_angle = atan2(eigvects[2, 1], eigvects[1, 1]) * (180 / pi());
}