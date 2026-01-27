library(terra)

## UTM
r <- rast(ext(c(1e+05, 9e+05, 0, 9e+06)), res = 1e5, crs = "EPSG:26918")
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))

maps::map(add = TRUE)


## OMERC
r <- rast(
  ext(c(-1e6, 1e06, -1.5e6, 3e6)),
  res = 1e5,
  crs = "+proj=omerc +lonc=-75.5 +lat_0=35.2 +gamma=-35.02"
)
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))

maps::map(add = TRUE)

r <- rast(
  ext(c(-5e6, 5e06, -1.5e6, 3e6)),
  res = 1e5,
  crs = "+proj=omerc +lonc=-75.5 +lat_0=35.2 +gamma=-35.02"
)
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))

maps::map(add = TRUE)


r <- rast(
  ext(c(-5e6, 5e06, -1.5e6, 3e6)),
  res = 1e5,
  crs = "+proj=omerc +lonc=-75.5 +lat_0=35.2 +gamma=-35.02"
)
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))

maps::map(add = TRUE)

## MERC
r <- rast(
  ext(c(-8e6, -5.2e06, 2.3e6, 5e6)),
  res = 1e5,
  crs = "+proj=merc +lat_ts=35.2"
)
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))
maps::map(add = TRUE)


r <- rast(ext(c(-9.2e6, -7e06, 3e6, 7.6e6)), res = 1e5, crs = "+proj=merc")
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))
maps::map(add = TRUE)


## AEQD

r <- rast(
  ext(c(-1e6, 1e6, -1e6, 1.5e6)),
  res = 1e5,
  crs = "+proj=aeqd +lat_0=35.2 +lon_0=-75.5"
)
plot(project(densify(as.polygons(r), 5000), "EPSG:4326"))
maps::map(add = TRUE)
