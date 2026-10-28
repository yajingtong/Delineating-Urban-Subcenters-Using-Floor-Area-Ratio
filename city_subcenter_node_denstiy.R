#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(sf); library(dplyr); library(readr)})
road_classes <- c("motorway", "trunk", "primary", "secondary", "tertiary", "unclassified",
  "residential", "living_street", "service", "motorway_link", "trunk_link", "primary_link",
  "secondary_link", "tertiary_link")

build_road_graph <- function(roads) {
  stopifnot(!is.na(st_crs(roads)))
  roads <- st_transform(roads, 2154)
  roads <- roads[!st_is_empty(roads), ]
  if (!nrow(roads)) stop("Road layer is empty")
  if (!all(as.character(st_geometry_type(roads)) %in% c("LINESTRING", "MULTILINESTRING")))
    stop("Roads must contain only LINESTRING or MULTILINESTRING geometries")
  coords <- st_coordinates(suppressWarnings(st_cast(st_geometry(roads), "LINESTRING")))
  if (nrow(coords) < 2L || any(!is.finite(coords[, c("X", "Y")])))
    stop("Road coordinates are missing or invalid")
  keys <- paste(as.integer(round(coords[, "X"] * 100)),
                as.integer(round(coords[, "Y"] * 100)), sep = ":")
  first <- !duplicated(keys)
  node_xy <- coords[first, c("X", "Y"), drop = FALSE]
  ids <- match(keys, keys[first])
  i <- seq_len(nrow(coords) - 1L)
  i <- i[coords[i, "L1"] == coords[i + 1L, "L1"] & ids[i] != ids[i + 1L]]
  if (!length(i)) stop("Road layer has no nonzero edges")
  edges <- cbind(ids[i], ids[i + 1L])
  weights <- sqrt(rowSums((coords[i, c("X", "Y"), drop = FALSE] -
                            coords[i + 1L, c("X", "Y"), drop = FALSE])^2))
  stopifnot(all(is.finite(weights)), all(weights > 0), !anyNA(edges))
  graph <- igraph::make_empty_graph(n = nrow(node_xy), directed = FALSE)
  graph <- igraph::add_edges(graph, as.vector(t(edges)))
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                            edge.attr.comb = list(weight = "min", "ignore"))
  list(graph = graph, node_xy = node_xy)
}

local_node_measures <- function(patches, graph, node_xy) {
  patches <- st_transform(patches, 2154)
  if (!all(st_is_valid(patches)) || any(st_is_empty(patches)))
    stop("Subcentre polygons must be valid and nonempty")
  if (!all(as.character(st_geometry_type(patches)) %in% c("POLYGON", "MULTIPOLYGON")))
    stop("Subcentres must be polygons")
  required <- c("city", "patch_label", "hh_patch_id", "rcoi_type")
  if (!all(required %in% names(patches))) stop("Missing subcentre identity/type fields")
  if (anyNA(st_drop_geometry(patches)[, required[1:3]]) ||
      anyDuplicated(patches$hh_patch_id) || anyDuplicated(patches$patch_label))
    stop("City subcentre identifiers must be complete and unique")
  stopifnot(igraph::vcount(graph) == nrow(node_xy), !igraph::is_directed(graph))
  degree <- igraph::degree(graph, loops = FALSE)
  components <- igraph::components(graph)
  largest <- which.max(components$csize)
  if (!any(degree > 0)) stop("No active road nodes")
  candidate <- rep(FALSE, nrow(node_xy))
  for (i in seq_len(nrow(patches))) {
    bb <- st_bbox(patches[i, ])
    candidate <- candidate | (node_xy[, 1] >= bb["xmin"] & node_xy[, 1] <= bb["xmax"] &
                              node_xy[, 2] >= bb["ymin"] & node_xy[, 2] <= bb["ymax"])
  }
  ids <- which(candidate & degree > 0)
  nodes <- st_as_sf(data.frame(graph_node_id = ids, x = node_xy[ids, 1], y = node_xy[ids, 2],
    degree = degree[ids], component_id = components$membership[ids]), coords = c("x", "y"), crs = 2154)
  hits <- st_intersects(nodes, patches)
  # Preserve the Paris rule: never silently double-count or arbitrarily assign shared nodes.
  if (any(lengths(hits) > 1L)) stop("Road nodes intersect multiple subcentres; resolve polygon overlap/shared boundaries")
  member <- lengths(hits) == 1L
  nodes <- nodes[member, ]
  nodes$patch_index <- vapply(hits[member], function(x) x[1], integer(1))
  nodes$city <- patches$city[nodes$patch_index]
  nodes$subcentre <- patches$patch_label[nodes$patch_index]
  nodes$patch_id <- patches$hh_patch_id[nodes$patch_index]
  area <- as.numeric(st_area(patches)) / 1e6
  if (any(!is.finite(area) | area <= 0)) stop("Nonpositive polygon area")
  metrics <- bind_rows(lapply(seq_len(nrow(patches)), function(i) {
    z <- nodes[nodes$patch_index == i, ]
    deg <- z$degree
    street <- deg != 2; junction <- deg >= 3
    data.frame(city = patches$city[i], subcentre = patches$patch_label[i],
      patch_id = patches$hh_patch_id[i], functional_type = patches$rcoi_type[i], area_km2 = area[i],
      all_graph_vertices = length(deg), geometry_vertices_degree2 = sum(deg == 2),
      street_nodes_excluding_degree2 = sum(street), junctions_degree3plus = sum(junction),
      dead_ends_degree1 = sum(deg == 1), four_plus_way_junctions = sum(deg >= 4),
      street_node_density_per_km2 = sum(street) / area[i], junction_density_per_km2 = sum(junction) / area[i],
      mean_street_node_degree = if (any(street)) mean(deg[street]) else NA_real_,
      mean_junction_degree = if (any(junction)) mean(deg[junction]) else NA_real_,
      dead_end_fraction_of_street_nodes = if (any(street)) sum(deg == 1) / sum(street) else NA_real_,
      global_components_present = length(unique(z$component_id)),
      share_graph_vertices_in_largest_component = if (length(deg)) mean(z$component_id == largest) else NA_real_)
  }))
  stopifnot(sum(metrics$all_graph_vertices) == nrow(nodes))
  summary <- data.frame(city = patches$city[1], subcentres = nrow(patches),
    graph_vertices = igraph::vcount(graph), graph_edges = igraph::ecount(graph),
    nonisolated_vertices = sum(degree > 0), connected_components = components$no,
    largest_component_vertices = max(components$csize),
    largest_component_fraction = max(components$csize) / sum(degree > 0),
    subcentres_without_road_nodes = sum(metrics$all_graph_vertices == 0))
  list(metrics = metrics, nodes = nodes, summary = summary, patches = patches)
}

# Cache signature intentionally matches the original Paris graph implementation.
load_city_graph <- function(path, layer, cache) {
  path <- normalizePath(path, mustWork = TRUE)
  signature <- list(path = path, size = file.info(path)$size,
    mtime = as.numeric(file.info(path)$mtime), crs = 2154, rounding_m = .01, version = 1)
  if (file.exists(cache)) {
    old <- readRDS(cache)
    # New caches also record the layer; the existing Paris cache predates that field.
    legacy_paris <- is.null(old$layer) && basename(path) == "paris_osm_roads_zipf_boundary.gpkg" &&
      layer == "osm_roads_paris_zipf_boundary"
    if (identical(old$signature, signature) && (identical(old$layer, layer) || legacy_paris)) return(old)
    rm(old); gc(FALSE)
  }
  roads <- st_read(path, layer = layer, quiet = TRUE)
  if (!"highway" %in% names(roads)) stop("Road layer requires a highway field to check road-class comparability")
  types <- as.character(roads$highway)
  if (anyNA(types) || any(!types %in% road_classes))
    stop("Road layer contains missing/extra highway classes; prepare the same 14 Paris road classes first")
  result <- build_road_graph(roads)
  result$signature <- signature; result$layer <- layer; result$source_features <- nrow(roads)
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, cache, compress = TRUE)
  result
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  allowed <- c("root", "patches", "manifest", "output", "cities", "prepare")
  if (any(!grepl("^--[a-z]+(=|$)", args))) stop("Use --key=value arguments; see README.md")
  keys <- sub("=.*$", "", sub("^--", "", args))
  if (any(!keys %in% allowed) || anyDuplicated(keys)) stop("Unknown or repeated argument")
  opts <- setNames(lapply(args, function(x) if (grepl("=", x, fixed = TRUE)) sub("^[^=]*=", "", x) else TRUE), keys)
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  default_root <- if (length(file_arg)) file.path(dirname(sub("^--file=", "", file_arg[1])), "../..") else getwd()
  root <- normalizePath(if (is.null(opts$root)) default_root else opts$root, mustWork = TRUE)
  resolve <- function(p) if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) p else file.path(root, p)
  value <- function(key, default) if (is.null(opts[[key]])) default else as.character(opts[[key]])
  output <- resolve(value("output", "outputs/all_cities_gt5_local_node_connectivity"))
  patch_path <- resolve(value("patches", "outputs/subcentre_standardized_poi_rcoi_32_cities/All_32_cities_standardized_POI_RCOI.gpkg"))
  manifest_path <- resolve(value("manifest", "scripts/all_cities_local_node_connectivity/road_sources.csv"))
  for (pkg in c("sf", "dplyr", "readr", "igraph", "ggplot2"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the required R package: ", pkg)
  patches <- st_read(patch_path, layer = "classified_patches", quiet = TRUE)
  eligible <- st_drop_geometry(patches) |> count(city, name = "subcentres") |>
    filter(subcentres > 5) |> arrange(desc(subcentres), city)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  write_csv(eligible, file.path(output, "eligible_cities.csv"))
  if (!is.null(opts$cities)) {
    requested <- strsplit(as.character(opts$cities), ",", fixed = TRUE)[[1]]
    if (any(!requested %in% eligible$city)) stop("--cities includes unknown cities or cities with <=5 subcentres")
    eligible <- filter(eligible, city %in% requested)
  }
  manifest <- read_csv(manifest_path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  required <- c("city", "road_path", "road_layer", "graph_cache_path")
  if (!all(required %in% names(manifest)) || anyNA(manifest$city) || anyDuplicated(manifest$city))
    stop("Manifest requires unique city rows and columns: ", paste(required, collapse = ", "))
  if (isTRUE(opts$prepare)) {
    message("Eligible cities: ", paste(eligible$city, collapse = ", "), ". Road inputs: ", manifest_path)
    return(invisible(eligible))
  }
  results <- list(); summaries <- list(); statuses <- list()
  for (city in eligible$city) {
    message("Processing ", city)
    row <- manifest[manifest$city == city, ]
    status <- tryCatch({
      if (nrow(row) != 1 || anyNA(row[, required])) stop("Missing road manifest entry or required field")
      road_path <- resolve(row$road_path)
      if (!file.exists(road_path) || file.info(road_path)$size <= 0) stop("Road source unavailable: ", road_path)
      graph <- load_city_graph(road_path, row$road_layer, resolve(row$graph_cache_path))
      city_patches <- patches[patches$city == city, ]
      city_patches <- city_patches[order(as.integer(sub("^C", "", city_patches$patch_label))), ]
      z <- local_node_measures(city_patches, graph$graph, graph$node_xy)
      city_out <- file.path(output, city); dir.create(city_out, recursive = TRUE, showWarnings = FALSE)
      z$summary$source_road_features <- graph$source_features
      z$summary$road_path <- normalizePath(road_path); z$summary$road_layer <- row$road_layer
      write_csv(z$metrics, file.path(city_out, "subcentre_local_node_measures.csv"))
      write_csv(z$summary, file.path(city_out, "network_summary.csv"))
      gpkg <- file.path(city_out, "subcentre_local_node_connectivity.gpkg")
      spatial <- left_join(z$patches, z$metrics, by = c("city", "patch_label" = "subcentre", "hh_patch_id" = "patch_id"))
      st_write(spatial, gpkg, layer = "subcentre_connectivity", delete_layer = TRUE, quiet = TRUE)
      if (nrow(z$nodes)) st_write(z$nodes, gpkg, layer = "network_vertices_inside_subcentres", delete_layer = TRUE, quiet = TRUE)
      p <- ggplot2::ggplot(z$metrics, ggplot2::aes(x = reorder(subcentre, junction_density_per_km2), y = junction_density_per_km2)) +
        ggplot2::geom_col(fill = "#B2182B") + ggplot2::coord_flip() +
        ggplot2::labs(title = city, x = "Subcentre", y = "Junctions (degree >= 3) per km²") +
        ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
      ggplot2::ggsave(file.path(city_out, "junction_density.png"), p, width = 20, height = 28, units = "cm", dpi = 180, bg = "white")
      ggplot2::ggsave(file.path(city_out, "junction_density.pdf"), p, width = 20, height = 28, units = "cm", useDingbats = FALSE)
      results[[city]] <- z$metrics; summaries[[city]] <- z$summary
      rm(graph, z); gc(FALSE)
      data.frame(city = city, status = "complete", detail = city_out)
    }, error = function(e) {
      message(city, ": ", conditionMessage(e))
      data.frame(city = city, status = "unavailable_or_failed", detail = conditionMessage(e))
    })
    statuses[[city]] <- status
    write_csv(bind_rows(statuses), file.path(output, "run_status.csv"))
  }
  # These aggregates represent only this invocation; status lists every requested city.
  if (length(results)) {
    write_csv(bind_rows(results), file.path(output, "all_completed_cities_local_node_measures.csv"))
    write_csv(bind_rows(summaries), file.path(output, "all_completed_cities_network_summary.csv"))
  } else {
    for (name in c("all_completed_cities_local_node_measures.csv", "all_completed_cities_network_summary.csv"))
      if (file.exists(file.path(output, name))) unlink(file.path(output, name))
  }
  message("Saved results: ", output)
  if (any(bind_rows(statuses)$status != "complete"))
    warning("Some eligible cities lack usable inputs or failed validation; inspect run_status.csv", call. = FALSE)
  invisible(bind_rows(statuses))
}
if (sys.nframe() == 0L) main()
