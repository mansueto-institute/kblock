| Variable | Description | Type | Unit | 
|---|---|---|---|
| building | Type of building | str |  |
| highway | Method of accessing building | str |  |
| name | Building name, if there is one | str |  |
| natural | Whether the feature (in this case, a building), is man-made or not. Should always be yes. | str |  |
| osm_id | Building OSM ID | int |  |
| osm_way_id | Way OSM ID (a Way in OSM can be a line or an area i.e. a park) | int |  |
| other_tags | Various data about the building | Comma separated str |  |
| waterway | Type of waterway, if any. Should be None. | str |  |
| gadm_code | Most local GADM the building and block belong to | str |  |
| bldg_id | ID # of the building within the block | str |  |
| bldg_pop | Estimated population count of the building | float | People |
| block_area | The area of the block that the building belongs to | float | km^2 |
| block_bldg_count | The number of buildings in the block | float | Building |
| block_bldg_area  | The summed areas of the block's buildings | float | km^2 |
| block_bldg_area_density | The ratio of block_bldg_area to block_area | float | |
| block_bldg_count_density | The ratio of block_bldg_count to block_area | float| Buildings / km^2 | 
| block_pop | The summed population estimate of the block | float | People |
| block_pop_density | The ratio of block_pop to block_area | float | People / km^2 |
| geometry | block geometry | shapely.Polygon |  | 
