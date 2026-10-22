#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(sf); library(dplyr); library(readr); library(ggplot2)})

arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(arg)) normalizePath(file.path(dirname(sub("^--file=", "", arg[1])), "../..")) else getwd()
out <- file.path(root, "outputs/paris_subcentre_network_connectivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
road_path <- file.path(root, "outputs/osm_roads_paris_zipf_boundary/paris_osm_roads_zipf_boundary.gpkg")
patch_path <- file.path(root, "outputs/subcentre_standardized_poi_rcoi_32_cities/by_city/Paris_standardized_POI_RCOI.gpkg")
context <- file.path(root, "outputs/all_41_cities_buildings_subcentre_types_updated/all_41_cities_updated_functional_subcentres.gpkg")
patches <- st_read(patch_path, layer = "classified_patches", quiet = TRUE) |> st_transform(2154) |>
  mutate(patch_order = as.integer(sub("^C", "", patch_label))) |> arrange(patch_order)
center <- st_read(context, query = "SELECT * FROM urban_centres WHERE city = 'Paris'", quiet = TRUE) |> st_transform(2154)
stopifnot(nrow(patches) == 44, nrow(center) == 1, all(st_is_valid(patches)))
centroids <- suppressWarnings(st_centroid(patches))
cache <- file.path(out, "road_graph_cache.rds")
signature <- list(path = road_path, size = file.info(road_path)$size,
                  mtime = as.numeric(file.info(road_path)$mtime), crs = 2154, rounding_m = .01, version = 1)
cached <- if (file.exists(cache)) readRDS(cache) else NULL
if (!is.null(cached) && identical(cached$signature, signature)) {
  graph <- cached$graph; node_xy <- cached$node_xy; source_features <- cached$source_features
  message("Reusing checked road-graph cache")
  rm(cached)
} else {
  rm(cached)

roads <- st_read(road_path, query = "SELECT geom FROM osm_roads_paris_zipf_boundary", quiet = TRUE) |> st_transform(2154)
  source_features <- nrow(roads)
  stopifnot(source_features == 661838)
  lines <- st_cast(st_geometry(roads), "LINESTRING")
  coords <- st_coordinates(lines)
  rm(roads, lines); gc(FALSE)
  message("Building a shared-vertex graph from ", format(nrow(coords), big.mark = ","), " road vertices")
  keys <- paste(as.integer(round(coords[, "X"] * 100)), as.integer(round(coords[, "Y"] * 100)), sep = ":")
  first <- !duplicated(keys)
  node_xy <- coords[first, c("X", "Y"), drop = FALSE]
  ids <- match(keys, keys[first])
  i <- seq_len(nrow(coords) - 1)
  keep <- coords[i, "L1"] == coords[i + 1, "L1"] & ids[i] != ids[i + 1]
  i <- i[keep]
  edge_matrix <- cbind(ids[i], ids[i + 1])
  weights <- sqrt(rowSums((coords[i, c("X", "Y"), drop = FALSE] - coords[i + 1, c("X", "Y"), drop = FALSE])^2))
  stopifnot(all(weights > 0), !anyNA(edge_matrix))
  graph <- igraph::make_empty_graph(n = nrow(node_xy), directed = FALSE)
  graph <- igraph::add_edges(graph, as.vector(t(edge_matrix)))
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                            edge.attr.comb = list(weight = "min", "ignore"))
  rm(coords, keys, first, ids, i, keep, edge_matrix, weights); gc(FALSE)
  message("Saving reusable road graph: ", igraph::vcount(graph), " vertices; ", igraph::ecount(graph), " edges")
  saveRDS(list(signature = signature, graph = graph, node_xy = node_xy, source_features = source_features), cache, compress = TRUE)
}
degree <- igraph::degree(graph, loops = FALSE)
components <- igraph::components(graph)
largest <- which.max(components$csize)
active_ids <- which(degree > 0)
stopifnot(length(active_ids) > 0)
message("Assigning road nodes to the 44 subcentre polygons")
# Avoid constructing millions of sf points: first retain nodes in any patch bbox.
candidate <- rep(FALSE, nrow(node_xy))
for (i in seq_len(nrow(patches))) {
  bb <- st_bbox(patches[i, ])
  candidate <- candidate | (node_xy[, 1] >= bb["xmin"] & node_xy[, 1] <= bb["xmax"] &
                            node_xy[, 2] >= bb["ymin"] & node_xy[, 2] <= bb["ymax"])
}
candidate_ids <- which(candidate & degree > 0)
node_points <- st_as_sf(data.frame(graph_node_id = candidate_ids, x = node_xy[candidate_ids, 1],
  y = node_xy[candidate_ids, 2], degree = degree[candidate_ids],
  component_id = components$membership[candidate_ids]), coords = c("x", "y"), crs = 2154)
hits <- st_intersects(node_points, patches)
stopifnot(!any(lengths(hits) > 1))
member <- lengths(hits) == 1
node_points <- node_points[member, ]
node_points$patch_index <- vapply(hits[member], function(x) x[1], integer(1))
node_points$subcentre <- patches$patch_label[node_points$patch_index]
node_points$patch_id <- patches$hh_patch_id[node_points$patch_index]
area_km2 <- as.numeric(st_area(patches)) / 1e6
metrics <- bind_rows(lapply(seq_len(nrow(patches)), function(i) {
  z <- node_points[node_points$patch_index == i, ]
  deg <- z$degree
  # Degree-2 geometry vertices are not street intersections.
  street_nodes <- deg != 2
  junctions <- deg >= 3
  data.frame(subcentre = patches$patch_label[i], patch_id = patches$hh_patch_id[i],
    functional_type = patches$rcoi_type[i], area_km2 = area_km2[i],
    all_graph_vertices = length(deg), geometry_vertices_degree2 = sum(deg == 2),
    street_nodes_excluding_degree2 = sum(street_nodes), junctions_degree3plus = sum(junctions),
    dead_ends_degree1 = sum(deg == 1), four_plus_way_junctions = sum(deg >= 4),
    street_node_density_per_km2 = sum(street_nodes) / area_km2[i],
    junction_density_per_km2 = sum(junctions) / area_km2[i],
    mean_street_node_degree = if (any(street_nodes)) mean(deg[street_nodes]) else NA_real_,
    mean_junction_degree = if (any(junctions)) mean(deg[junctions]) else NA_real_,
    dead_end_fraction_of_street_nodes = if (any(street_nodes)) sum(deg == 1) / sum(street_nodes) else NA_real_,
    global_components_present = length(unique(z$component_id)),
    share_graph_vertices_in_largest_component = if (length(deg)) mean(z$component_id == largest) else NA_real_)
}))
stopifnot(nrow(metrics) == 44, sum(metrics$all_graph_vertices) == nrow(node_points))
message("Snapping 44 polygon centroids and the CBD to their nearest active road nodes")
anchor_xy <- rbind(st_coordinates(center), st_coordinates(centroids))
anchor_labels <- c("CBD", patches$patch_label)
nearest <- FNN::get.knnx(node_xy[active_ids, , drop = FALSE], anchor_xy, k = 1)
anchor_nodes <- active_ids[nearest$nn.index[, 1]]
snapping_m <- nearest$nn.dist[, 1]
anchor_table <- data.frame(location = anchor_labels, graph_node_id = anchor_nodes,
  centroid_easting_m = anchor_xy[, 1], centroid_northing_m = anchor_xy[, 2],
  snapped_node_easting_m = node_xy[anchor_nodes, 1], snapped_node_northing_m = node_xy[anchor_nodes, 2],
  snap_distance_m = snapping_m, node_degree = degree[anchor_nodes],
  component_id = components$membership[anchor_nodes],
  in_largest_component = components$membership[anchor_nodes] == largest)
write_csv(anchor_table, file.path(out, "centroid_and_CBD_node_snapping.csv"))
message("Computing weighted shortest paths between 45 snapped locations")
network_m <- igraph::distances(graph, v = anchor_nodes, to = anchor_nodes,
  mode = "all", weights = igraph::E(graph)$weight, algorithm = "dijkstra")
euclidean_km <- as.matrix(dist(anchor_xy)) / 1000
stopifnot(max(abs(network_m[is.finite(network_m)] - t(network_m)[is.finite(network_m)])) < 1e-5,
          all(diag(network_m) == 0))
dimnames(network_m) <- list(anchor_labels, anchor_labels)
matrix_km <- network_m / 1000; matrix_km[!is.finite(matrix_km)] <- NA_real_
write.csv(matrix_km, file.path(out, "road_distance_matrix_between_snapped_nodes_km.csv"), na = "NA")
ij <- which(upper.tri(network_m[-1, -1]), arr.ind = TRUE) + 1
pair_table <- data.frame(subcentre_1 = anchor_labels[ij[, 1]], subcentre_2 = anchor_labels[ij[, 2]],
  connected = is.finite(network_m[ij]), network_distance_between_snapped_nodes_km = matrix_km[ij],
  centroid_euclidean_distance_km = euclidean_km[ij],
  centroid_1_snap_distance_m = snapping_m[ij[, 1]], centroid_2_snap_distance_m = snapping_m[ij[, 2]],
  distance_with_straight_snap_connectors_km = matrix_km[ij] + (snapping_m[ij[, 1]] + snapping_m[ij[, 2]]) / 1000)
stopifnot(nrow(pair_table) == choose(44, 2))
metrics$centroid_snap_distance_m <- snapping_m[-1]
metrics$snapped_node_component <- anchor_table$component_id[-1]
metrics$snapped_node_in_largest_component <- anchor_table$in_largest_component[-1]
metrics$network_connected_to_CBD <- is.finite(network_m[-1, 1])
metrics$network_distance_to_CBD_km <- matrix_km[-1, 1]
metrics$centroid_euclidean_distance_to_CBD_km <- euclidean_km[-1, 1]
metrics$CBD_snap_distance_m <- snapping_m[1]
metrics$distance_to_CBD_with_straight_snap_connectors_km <- matrix_km[-1, 1] + (snapping_m[-1] + snapping_m[1]) / 1000
metrics$reachable_other_subcentres <- rowSums(is.finite(network_m[-1, -1])) - 1
metrics$mean_network_distance_to_reachable_subcentres_km <- vapply(seq_len(44), function(i) {
  x <- matrix_km[-1, -1][i, -i]; if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
}, numeric(1))
write_csv(metrics, file.path(out, "paris_44_subcentres_network_connectivity.csv"))
write_csv(pair_table, file.path(out, "paris_946_subcentre_network_distances_km.csv"))
write_csv(data.frame(source_road_features = source_features, graph_vertices = igraph::vcount(graph),
  graph_edges = igraph::ecount(graph), nonisolated_vertices = length(active_ids),
  connected_components = components$no, largest_component_vertices = max(components$csize),
  largest_component_fraction = max(components$csize) / length(active_ids),
  subcentres = 44, subcentres_connected_to_CBD = sum(metrics$network_connected_to_CBD),
  connected_subcentre_pairs = sum(pair_table$connected), total_subcentre_pairs = nrow(pair_table)),
  file.path(out, "network_summary.csv"))
patch_metrics <- left_join(patches, metrics, by = c("patch_label" = "subcentre", "hh_patch_id" = "patch_id"))
gpkg <- file.path(out, "paris_subcentre_network_connectivity.gpkg")
st_write(patch_metrics, gpkg, layer = "subcentre_connectivity", delete_layer = TRUE, quiet = TRUE)
st_write(node_points, gpkg, layer = "network_vertices_inside_subcentres", delete_layer = TRUE, quiet = TRUE)
anchors <- st_as_sf(anchor_table, coords = c("snapped_node_easting_m", "snapped_node_northing_m"), crs = 2154, remove = FALSE)
st_write(anchors, gpkg, layer = "snapped_centroids_and_CBD", delete_layer = TRUE, quiet = TRUE)
p <- ggplot(metrics, aes(x = reorder(subcentre, junction_density_per_km2), y = junction_density_per_km2)) +
  geom_col(fill = "#B2182B") + coord_flip() +
  labs(title = "Paris subcentre street-network connectivity", x = "Subcentre", y = "Junctions (degree >= 3) per km²",
    caption = "Degree is measured on the full undirected road graph before assigning nodes to subcentre polygons.\nDegree-2 geometry vertices are excluded from junction counts. Existing road geometries determine topology.") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0))
ggsave(file.path(out, "paris_subcentre_junction_density.png"), p, width = 20, height = 28, units = "cm", dpi = 180, bg = "white")
ggsave(file.path(out, "paris_subcentre_junction_density.pdf"), p, width = 20, height = 28, units = "cm",
  device = grDevices::pdf, encoding = "WinAnsi.enc", useDingbats = FALSE, bg = "white")
writeLines(c(
  "Paris node-based road connectivity: 44 existing standardized-POI subcentres.",
  paste("Road source:", road_path), paste("Patch source:", patch_path),
  "Use the complete 99.56 km Paris study-area road dataset, not the old Paris-commune node file.",
  "Road classes are the source's 14 highway categories, including service roads and motorway links; walking paths and railways are not added.",
  "Build an undirected graph by connecting consecutive existing LINESTRING vertices. Shared coordinates are quantized to a 1 cm grid in Lambert-93 (EPSG:2154).",
  "No new nodes are inserted at geometric crossings: shared existing geometry vertices define connections. Missing OSM node IDs mean topology and grade separation cannot be guaranteed; exact coincident geometry vertices may be conflated.",
  "One-way restrictions, access restrictions, barriers, turn restrictions and travel speeds are not modeled. Results describe 2D geometric network connectivity, not validated driving or pedestrian routes.",
  "Parallel duplicate node-to-node edges are collapsed to minimum length; loops are removed. Lengths use the original metric segment coordinates.",
  "Degree-2 shape vertices remain in the routing graph but are excluded from street-node and junction metrics. Junctions have degree >=3; dead ends have degree 1.",
  "Degrees and components are computed on the whole road graph, then nodes are assigned by polygon intersection; roads are not first clipped to each patch.",
  "Polygon areas, not complete-city areas or building-buffer areas, normalize node densities. Junction density is the primary local connectivity indicator.",
  "Subcentre polygons and CBD are transformed to EPSG:2154; true polygon centroids are snapped to the nearest nonisolated graph vertex using Euclidean nearest-neighbor search.",
  "Network distances are length-weighted shortest paths between snapped vertices. Off-network snap offsets are reported separately.",
  "The columns including straight snap connectors add those offsets for comparison; those straight connectors are not surveyed or routable road segments.",
  "Disconnected locations have NA network distances and explicit connected flags; no forced reconnection or largest-component fallback is applied.",
  "946 unique unordered subcentre pairs plus 44 CBD connections are exported. A 45-location distance matrix includes the CBD.",
  "Road data are clipped to the study boundary, so routes outside that extent cannot be evaluated. Source retrieval metadata are preserved in outputs/osm_roads_paris_zipf_boundary/.",
  "Node counts depend on source geometry quality and junction representation; nearby roundabout or dual-carriageway junctions are not consolidated."
), file.path(out, "README.txt"))
message("COMPLETE: local node metrics for 44 subcentres; ", sum(pair_table$connected), "/946 subcentre pairs connected; ",
        sum(metrics$network_connected_to_CBD), "/44 subcentres connected to CBD. Outputs: ", out)
