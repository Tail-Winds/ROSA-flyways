#!/usr/bin/env python3
# /// script
# dependencies = [
#   "pandas",
#   "xarray",
#   "zarr",
#   "dask",
#   "s3fs",
#   "fsspec",
# ]
# ///

import sys
import argparse
import pandas as pd
import xarray as xr
import dask

def main():
    parser = argparse.ArgumentParser(description="Extract JPL MUR SST v4.1 data for spatial points.")
    parser.add_argument("-i", "--input", default=sys.stdin, help="Input CSV file path (default: stdin)")
    parser.add_argument("-o", "--output", default=sys.stdout, help="Output CSV file path (default: stdout)")
    parser.add_argument("-w", "--workers", type=int, default=16, help="Number of Dask thread workers (default: 16)")
    args = parser.parse_args()

    # 1. Read input dataset
    df = pd.read_csv(args.input)
    df['date'] = pd.to_datetime(df["date"])

    # 2. Open S3 Zarr store
    sys.stderr.write("Opening S3 Zarr store...\n")
    ds = xr.open_zarr(
      "s3://mur-sst/zarr-v1",
      storage_options = {
        "anon": True,
        "client_kwargs": {"region_name": "us-west-2"},
        "default_fill_cache": False
      }, 
      consolidated = True
    )

    # Spatial bounding box (+/- 0.02 deg buffer)
    min_lat, max_lat = df['med_lat'].min() - 0.02, df['med_lat'].max() + 0.02
    min_lon, max_lon = df['med_lon'].min() - 0.02, df['med_lon'].max() + 0.02

    results = []

    for yr, sub_df in df.groupby(df['date'].dt.year):
        sub_df = sub_df.copy()
        unique_dates = sub_df['date'].drop_duplicates().sort_values().values
        
        sys.stderr.write(f"Processing {yr}: {len(unique_dates)} dates, {len(sub_df)} rows...\n")

        with dask.config.set(scheduler='threads', num_workers=args.workers):
            cube = (
                ds['analysed_sst']
                .sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
                .sel(time=unique_dates, method="nearest")
                .load()
            )

        extracted = cube.sel(
            time = xr.DataArray(sub_df['date'].values, dims="points"),
            lat = xr.DataArray(sub_df['med_lat'].values, dims="points"),
            lon = xr.DataArray(sub_df['med_lon'].values, dims="points"),
            method = 'nearest'
        )

        sub_df['sst_celsius'] = (extracted - 273.15).values
        results.append(sub_df)

    final_df = pd.concat(results, ignore_index=True)

    final_df['year'] = final_df['date'].dt.year
    final_df['wk'] = ((final_df['date'].dt.day_of_year - 1) // 7) + 1
    final_df = final_df.groupby(['med_lat', 'med_lon', 'year', 'wk']).agg(
        sst_celsius=('sst_celsius', 'mean')
    ).reset_index()


    final_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()