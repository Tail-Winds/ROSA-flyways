functions {
  // Equations from https://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
  real calc_k_prime(real lat, real lon, real lat0, real lon0) {
    // Eq 4 at the link above
    real cos_c = sin(lat0) * sin(lat) + cos(lat0) * cos(lat) * cos(lon - lon0);
    // Eq 3 at the link above
    real c = acos(fmin(1.0, fmax(-1.0, cos_c)));

    // Sets k' to 1 if c is very small to prevent dividing by 0
    return (c < 1e-9) ? 1.0 : c / sin(c);
  }

  matrix project_aeqd(vector lat, vector lon, real lat0, real lon0) {
    int N = rows(lat);
    matrix[N, 3] coords;

    for (n in 1:N) {
      real k_prime = calc_k_prime(lat[n], lon[n], lat0, lon0);
      // Eq. 1 (easting/X)
      coords[n, 1] = k_prime * cos(lat[n]) * sin(lon[n] - lon0);
      // Eq. 2 (norhing/Y)
      coords[n, 2] = k_prime * (cos(lat0) * sin(lat[n]) - 
        sin(lat0) * cos(lat[n]) * cos(lon[n] - lon0));

      coords[n, 3] = k_prime;
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
  real<lower=0, upper=10> sigma; 
}

transformed parameters {
  real lat0 = asin(center_vec[3]);
  real lon0 = atan2(center_vec[2], center_vec[1]);
}

model {
  for (n in 1:N) {
    // 1. Calculate distance c
    real cos_c = center_vec[1] * cos(lat_rad[n]) * cos(lon_rad[n]) + 
                 center_vec[2] * cos(lat_rad[n]) * sin(lon_rad[n]) + 
                 center_vec[3] * sin(lat_rad[n]);
    real c = acos(fmin(1.0, fmax(-1.0, cos_c)));
    
    // 2. Calculate k_prime (radial scale factor)
    real k_prime = (c < 1e-9) ? 1.0 : c / sin(c);
    
    // 3. Minimize |k_prime - 1| using Laplace Likelihood
    target += double_exponential_lpdf(k_prime | 1.0, sigma);
  }
  sigma ~ exponential(1);
}

generated quantities {
  matrix[N, 3] projected_points = project_aeqd(lat_rad, lon_rad, lat0, lon0);

  // Return median k_prime
  real k_prime = quantile(projected_points[, 3], 0.5);

  // --- PCA Logic ---
  matrix[2, 2] cov_mat;
  {
     matrix[N, 2] centered;
     for(j in 1:2) {
       centered[, j] = projected_points[, j] - mean(projected_points[, j]);
     }
     cov_mat = (centered' * centered) / (N - 1);
  }
  
  // Extract Eigenvectors (Principal Components)
  complex_matrix[2, 2] evecs_complex = eigenvectors(cov_mat);

  // 2. Extract the real part (since cov_matrix is symmetric, imaginary part is 0)
  matrix[2, 2] evecs = get_real(evecs_complex);

  // Angle of the first Principal Component (PC1)
  // This is the rotation needed to align the data's long axis
  real pca_angle_deg = atan2(evecs[2, 1], evecs[1, 1]) * (180 / pi());
}
