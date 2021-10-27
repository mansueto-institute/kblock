from pyarrow.parquet import ParquetFile
import geopandas as gpd
from typing import Dict, Optional, Sequence, Union
import json
from pathlib import Path

# For full reads and writes, use geopandas.read_parquet() and GeoDataFrame.to_parquet().
# geopandas.read_parquet() is also good for columnar subsets. 

# This file would need to be changed if we want to add it to geopandas. For our purposes, we don't care if 
# we asked for a geometry column or not, we'll just pass back whatever was requested. This function takes care of 
# that by passing back a dataframe if the geometry column is not in the set of columns. 
def read_row_groups(file, 
					row_groups: Sequence[str], 
					columns: Sequence[str]=None, 
					geometry: str='geometry', 
					use_threads: bool=True, 
					use_pandas_metadata: bool=False):
	pq_file = ParquetFile(file)
	pq_row_groups = pq_file.read_row_groups(row_groups=row_groups, columns=columns, use_threads=use_threads)


	if geometry in columns:
		# Copied code from https://github.com/geopandas/geopandas/blob/9607d9c96620093dbd12c346c61cea47b62b4dca/geopandas/io/arrow.py#L289
		df = pq_row_groups.to_pandas()
		metadata = pq_row_groups.schema.metadata
		if b"geo" not in metadata:
			raise ValueError(
				"""Missing geo metadata in Parquet/Feather file.
				Use pandas.read_parquet/read_feather() instead."""
			)

		try:
			metadata = decode_metadata(metadata.get(b"geo", b""))

		except (TypeError, json.decoder.JSONDecodeError):
			raise ValueError("Missing or malformed geo metadata in Parquet/Feather file")

		_validate_metadata(metadata)

		# Find all geometry columns that were read from the file.  May
		# be a subset if 'columns' parameter is used.
		geometry_columns = df.columns.intersection(metadata["columns"])

		if not len(geometry_columns):
			raise ValueError(
				"""No geometry columns are included in the columns read from
				the Parquet/Feather file.  To read this file without geometry columns,
				use pandas.read_parquet/read_feather() instead."""
			)

		geometry = metadata["primary_column"]

		# Missing geometry likely indicates a subset of columns was read;
		# promote the first available geometry to the primary geometry.
		if len(geometry_columns) and geometry not in geometry_columns:
			geometry = geometry_columns[0]

			# if there are multiple non-primary geometry columns, raise a warning
			if len(geometry_columns) > 1:
				warnings.warn(
					"Multiple non-primary geometry columns read from Parquet/Feather "
					"file. The first column read was promoted to the primary geometry."
				)

		# Convert the WKB columns that are present back to geometry.
		for col in geometry_columns:
			df[col] = from_wkb(df[col].values, crs=metadata["columns"][col]["crs"])

		return gpd.GeoDataFrame(df, geometry=geometry)
	else:
		df = pq_row_groups.to_pandas()
		return df


def read_row_group(file, 
				   row_group: str, 
				   columns: Sequence[str]=None, 
				   geometry: str='geometry', 
				   use_threads: bool=True, 
				   use_pandas_metadata: bool=False):
	pq_file = ParquetFile(file)
	pq_row_group = pq_file.read_row_group(row_groups=row_group, columns=columns, use_threads=use_threads)

	if geometry in columns:
		# Copied code from https://github.com/geopandas/geopandas/blob/9607d9c96620093dbd12c346c61cea47b62b4dca/geopandas/io/arrow.py#L289
		df = pq_row_groups.to_pandas()
		metadata = pq_row_group.schema.metadata
		if b"geo" not in metadata:
			raise ValueError(
				"""Missing geo metadata in Parquet/Feather file.
				Use pandas.read_parquet/read_feather() instead."""
			)

		try:
			metadata = decode_metadata(metadata.get(b"geo", b""))

		except (TypeError, json.decoder.JSONDecodeError):
			raise ValueError("Missing or malformed geo metadata in Parquet/Feather file")

		_validate_metadata(metadata)

		# Find all geometry columns that were read from the file.  May
		# be a subset if 'columns' parameter is used.
		geometry_columns = df.columns.intersection(metadata["columns"])

		if not len(geometry_columns):
			raise ValueError(
				"""No geometry columns are included in the columns read from
				the Parquet/Feather file.  To read this file without geometry columns,
				use pandas.read_parquet/read_feather() instead."""
			)

		geometry = metadata["primary_column"]

		# Missing geometry likely indicates a subset of columns was read;
		# promote the first available geometry to the primary geometry.
		if len(geometry_columns) and geometry not in geometry_columns:
			geometry = geometry_columns[0]

			# if there are multiple non-primary geometry columns, raise a warning
			if len(geometry_columns) > 1:
				warnings.warn(
					"Multiple non-primary geometry columns read from Parquet/Feather "
					"file. The first column read was promoted to the primary geometry."
				)

		# Convert the WKB columns that are present back to geometry.
		for col in geometry_columns:
			df[col] = from_wkb(df[col].values, crs=metadata["columns"][col]["crs"])

		return gpd.GeoDataFrame(df, geometry=geometry)
	else:
		df = pq_row_group.to_pandas()
		return df


def append_to_parquet(gdf: GeoDataFrame,
					  filepath: Union[str, Path]):
	pass


def decode_metadata(metadata_str):
    """Decode a UTF-8 encoded JSON string to dict
    Parameters
    ----------
    metadata_str : string (UTF-8 encoded)
    Returns
    -------
    dict
    """
    if metadata_str is None:
        return None

    return json.loads(metadata_str.decode("utf-8"))


def _validate_metadata(metadata):
    """Validate geo metadata.
    Must not be empty, and must contain the structure specified above.
    Raises ValueError if metadata is not valid.
    Parameters
    ----------
    metadata : dict
    """

    if not metadata:
        raise ValueError("Missing or malformed geo metadata in Parquet/Feather file")

    required_keys = ("primary_column", "columns")
    for key in required_keys:
        if metadata.get(key, None) is None:
            raise ValueError(
                "'geo' metadata in Parquet/Feather file is missing required key: "
                "'{key}'".format(key=key)
            )

    if not isinstance(metadata["columns"], dict):
        raise ValueError("'columns' in 'geo' metadata must be a dict")

    # Validate that geometry columns have required metadata and values
    required_col_keys = ("crs", "encoding")
    for col, column_metadata in metadata["columns"].items():
        for key in required_col_keys:
            if key not in column_metadata:
                raise ValueError(
                    "'geo' metadata in Parquet/Feather file is missing required key "
                    "'{key}' for column '{col}'".format(key=key, col=col)
                )

        if column_metadata["encoding"] != "WKB":
            raise ValueError("Only WKB geometry encoding is supported")
