library(qs2)
library(data.table)

# need position per individual per day
# STRIPED BASS
sb <- qs2::qs_read("data/RAW/hud_detects.qs2") |>
  as.data.table() |>
  _[, let(
    year = year(date.local),
    week = week(date.local)
  )] |>
  _[,
    .(
      mean_lat = mean(lat, na.rm = T),
      mean_lon = mean(long, na.rm = T),
      species = 'striped bass'
    ),
    by = c('week', 'year', 'transmitter')
  ]

fwrite(sb, "data/sb_weekly.csv")


# SAND TIGER
st <- fread("data/RAW/sandtigers.csv") |>
  _[, let(
    DateandTime = as.POSIXct(DateandTime, tz = "UTC", format = "%m/%d/%Y %H:%M")
  )] |>
  _[,
    let(
      week = week(DateandTime),
      year = year(DateandTime),
      transmitter = Transmitter
    )
  ] |>
  _[,
    .(
      mean_lat = mean(Latitude, na.rm = T),
      mean_lon = mean(Longitude, na.rm = T),
      species = 'sand tiger'
    ),
    by = c('week', 'year', 'transmitter')
  ]

fwrite(st, "data/st_weekly.csv")
