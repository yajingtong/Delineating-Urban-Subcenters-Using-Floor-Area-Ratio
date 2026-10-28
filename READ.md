# Local road-node connectivity for cities with more than five subcentres

## Run

From the project root:

```sh
Rscript scripts/all_cities_local_node_connectivity/compute_local_node_connectivity.R
```

Required R packages: sf, dplyr, readr, igraph, ggplot2. Packages are checked, never automatically installed. This script computes **local node measures**, without pairwise routing or modeled flows.

The script automatically counts rows in the standardized `classified_patches` layer and selects **strictly more than five** subcentres. Current eligible cities: Paris (44), Lyon (14), Marseille (10), Montpellier (8), Toulouse (7), Nantes (6). Cities with exactly five are excluded. Eligibility is recalculated on every run, rather than hardcoded.

`road_sources.csv` maps each city to its road GeoPackage, layer, and graph-cache path. All six entries point to existing usable inputs. Marseille, Montpellier, Toulouse and Nantes were extracted from Geofabrik OSM regional snapshots by `scripts/four_cities_osm_roads/retrieve_four_city_roads.R`; their source dates and extraction summaries are in `outputs/osm_roads_four_cities_retrieval/`. Add manifest rows if additional cities become eligible. Paths may be absolute or relative to the project root. Each road layer must describe that city's full study extent (including the space between subcentres), with the same 14 Paris highway classes. The script checks line geometries, CRS and class values; source completeness and extent must be checked when preparing new inputs. It deliberately does not infer full coverage from a bounding box, download roads, use Paris-commune nodes, or substitute empty networks. Lyon uses the existing full-study-area bounding-box extraction; retain this source extent in interpreting its component metrics.

Missing inputs and per-city validation errors are recorded in `run_status.csv`; other cities continue. Check this file before treating a batch as complete. A warning is emitted for any incomplete city. Results are saved under `outputs/all_cities_gt5_local_node_connectivity/` by default.

Optional commands:

```sh
# List eligible cities without constructing graphs.
Rscript scripts/all_cities_local_node_connectivity/compute_local_node_connectivity.R --prepare

# Run a subset; the >5 rule still applies.
Rscript scripts/all_cities_local_node_connectivity/compute_local_node_connectivity.R --cities=Paris,Lyon

# Override project root, standardized patches, manifest, or output directory.
Rscript scripts/all_cities_local_node_connectivity/compute_local_node_connectivity.R --root=/absolute/project/path --manifest=/absolute/road_sources.csv --output=outputs/custom_local_nodes
```

`--patches=path.gpkg` replaces the default standardized GeoPackage; it must contain `classified_patches`, with `city`, `patch_label`, `hh_patch_id`, `rcoi_type`, and polygon geometry. Each row represents one subcentre. `--prepare` also saves `eligible_cities.csv`. Only long `--key=value` arguments and `--prepare` are supported.

## Same local measures as Paris

For every subcentre the output includes:

- Polygon area in km²; all active graph vertices; degree-2 geometry vertices.
- Street-node count (degree != 2), junction count (degree >= 3), dead ends (degree = 1), and four-or-more-way junctions.
- Street-node and junction densities per polygon km².
- Mean street-node degree and mean junction degree.
- Dead ends divided by street-node count.
- Number of full-network components represented and fraction of polygon vertices in the full network's largest component.

Degrees and components are computed on the **full supplied city road graph**, before assigning nodes by polygon intersection. Road lines are not clipped to each subcentre. Graphs are undirected, built from consecutive existing line vertices after projection to EPSG:2154; shared coordinates use a 1 cm grid. Duplicate parallel edges collapse to minimum length; loops are removed. Degree-2 vertices remain in the graph but do not count as street intersections. Lengths use original metric segment coordinates. No nodes are inserted at geometric crossings. Missing OSM node IDs, grade separation, road representation, and source extent affect topology. One-way, access and turn restrictions, travel speeds, and roundabout consolidation are not modeled. Results describe geometric connectivity, not validated travel or observed traffic.

The 14 permitted highway values are motorway, trunk, primary, secondary, tertiary, unclassified, residential, living_street, service, and motorway_link, trunk_link, primary_link, secondary_link, tertiary_link. Road layers with other/missing classes stop with an explanation, rather than silently mixing definitions. Uncached road layers must expose `highway`. The original Paris graph cache is reused when its source signature matches. New caches also check the layer name; delete a cache to force rebuilding. Cache signatures use source path, size, modification time, CRS, rounding, and method version. Each cache is explicitly scoped to its city's manifest row.

An empty subcentre has zero counts/densities and undefined (NA) means, fractions and largest-component share. Multiple polygon matches for one road node cause that city to stop, following the original Paris rule; the script does not invent an assignment or double-count nodes. Invalid polygons require correction at the source. Where two largest components tie in size, igraph's first maximum component is used, as in Paris.

## Outputs

Each completed city's folder contains `subcentre_local_node_measures.csv`, `network_summary.csv`, a GeoPackage with `subcentre_connectivity` and `network_vertices_inside_subcentres` layers, and equally sized 20 × 28 cm junction-density charts (PNG/PDF, red bars). Original functional-type codes are preserved; `Mixed` is the existing code for RCOI.

Batch outputs are `eligible_cities.csv`, `run_status.csv`, `all_completed_cities_local_node_measures.csv`, and `all_completed_cities_network_summary.csv`. Aggregates contain only cities completed during the current invocation; `--cities` narrows the invocation. Existing per-city outputs can remain from earlier runs: always consult the current status file. No combined aggregate is retained if no city completes. Input roads, patches and existing Paris outputs are not modified; graph caches may be refreshed when their source signature changes.

## Verification

From the project root:

```sh
Rscript scripts/all_cities_local_node_connectivity/tests/verify_local_measures.R
```

Checks include a known graph with a shape vertex and duplicate segment, whole-network degree before polygon assignment, empty polygons, disconnected geometric crossings, and rejection of ambiguous assignments. When the original Paris cache/results are present, every local-measure column is compared for all 44 Paris subcentres at tolerance 1e-12. Verification does not require downloading data.


## Geographic node-density maps

After computing the four cities, run:

```sh
Rscript scripts/all_cities_local_node_connectivity/map_city_node_densities.R
```

The default cities are Marseille, Montpellier, Toulouse and Nantes. City names can be supplied as separate positional arguments. Map inputs are the completed subcentre GeoPackages, road-source manifest, existing INSEE boundaries and urban-centre points, and the existing all-city building display rasters. Additional cities require a `study_boundary` layer in their road GeoPackage and a building display raster.

Maps show both junction density (degree >=3/km²) and street-node density (degree !=2/km²), with labels showing each subcentre and rounded density. A common red scale is used across cities for each measure, with identical square geographic extents that include all road-study and INSEE boundaries. The current four-city comparison uses 140 × 140 km panels. The INSEE boundaries can extend beyond the circular road study areas; densities describe subcentres within the road study area. Road networks are rasterized at approximately 50 m solely for display; this does not alter any graph, node metric or polygon density. All building footprints available in the pre-existing display rasters are retained.

The outputs are in `outputs/all_cities_gt5_local_node_connectivity/four_city_density_maps/`: two combined comparisons (PNG/PDF/preview), individual city maps for both measures, the shared map extents and colour limits, and a CSV with all 31 subcentres. Full computed graphs and node coordinates are in `outputs/all_cities_gt5_local_node_connectivity/graph_cache/`; per-city node GeoPackages and original metric CSVs remain in each city folder. Maps show local density and geographic network context; they do not estimate intercity or subcentre flows.

## Subcentre road-connectivity maps

To map every unordered subcentre pair along the computed road graph:

```sh
Rscript scripts/all_cities_local_node_connectivity/map_subcentre_road_connections.R
```

Defaults to Marseille, Montpellier, Toulouse and Nantes. Supply city names as separate arguments to change the list. The script reuses each city's full-resolution graph cache and snaps true polygon centroids to their nearest nonisolated road vertices. Each unordered pair contributes one equal-unit Dijkstra shortest path. Road colour and width show the number of selected pair paths sharing a road segment. These are modeled potential connections, not observed trips or traffic counts. Tied alternatives are not split. CBD routes and straight off-network snap connectors are excluded. Current graphs are undirected and do not apply access, turn or one-way restrictions. Every subcentre anchor is connected in the current four-city inputs. The script stops if an input contains disconnected anchors or duplicate snapped anchors, rather than silently omitting connections.

Each city's `all_pair_flow_maps/` folder contains individual PNG, preview and PDF maps, route and shared-edge CSVs, and a GeoPackage containing all selected complete shortest paths, shared road segments and snapped subcentre nodes. Nearest-node snapping and the distance matrix are saved in each city's parent folder. Node snap distances are retained for interpretation. The full INSEE urban boundary, available building footprints, subcentre functions, CBD and C labels are included. `Mixed` is displayed as RCOI. The current maps use a city-specific full-boundary view and a city-specific count scale; raw overlap counts should be interpreted relative to that city's number of pairs.

The batch summary is `outputs/all_cities_gt5_local_node_connectivity/road_connectivity_maps/four_city_route_summary.csv`: Marseille has 45 pairs, Montpellier 28, Toulouse 21 and Nantes 15 (109 total). Every route length is independently checked against a weighted distance matrix within 1e-6 km, and shared-edge counts are checked against total traversed edge counts. Route caches check source signature, snapped anchor IDs and graph edge count before reuse.
