library(tidyverse)
library(polmineR)
library(cwbtools)
library(future)
library(hash)
library(rjson)
library(tidygraph)
library(visNetwork)


# install GermaParl
cwbtools::CORPUS_install(doi = "10.5281/zenodo.3742113")

"GERMAPARL" %>% kwic("Test")

# +++
# if the test worked, you need to export the corpus as a .vrt file (see the description on Github)
# now, proceed in Python
# +++


# back in R, switch on polmineR's multiprocessing
NUM_WORKERS <- availableCores() - 2
options("polmineR.mc" = TRUE)
options("polmineR.cores" = NUM_WORKERS)

# enter your corpus name
CORPUS <- "GERMA_CLUSTERS"

# map your metadata categories to the respective s_attribute
CORPUS %>% s_attributes
META_S_ATT <-
  list(
    "1990" = "decade",
    "2000" = "decade",
    "2010" = "decade",
    "CDU" = "party",
    "SPD" = "party",
    "FDP" = "party",
    "GRUENE" = "party",
    "LINKE" = "party"
  )

# map your metadata categories to the respective p_attribute
CORPUS %>% p_attributes
META_P_ATT <-
  list(
    "1990" = "dcluster",
    "2000" = "dcluster",
    "2010" = "dcluster",
    "CDU" = "pcluster",
    "SPD" = "pcluster",
    "FDP" = "pcluster",
    "GRUENE" = "pcluster",
    "LINKE" = "pcluster"
  )

# if you want to use specific reference corpora, map your metadata categories to the reference corpus metadata here.
# (especially useful for diachronic analyses)
# if you enter a mapping, it is assumed that your subcorpus is not contained in the reference corpus
# don't enter anything if you want your subcorpora to be compared to the whole corpus
KEYNESS_COMPARISONS <-
  list(
    "2000" = "1990",
    "2010" = "2000"
    )

# generated in Python
META_CLUSTER_WORDS <-
  fromJSON(file = "path/output/meta_cluster_words.json")

# minimum cosine similarity to thematic centroid to be included as a node in the graph
T_SIM <- 0.25

#statistical measure for calculation of keywords
METHOD_KEYWORDS <- "ll"

# statistical measure for calculation of collocations
METHOD_COLLOCS <- "ll"

# test your enriched corpus
test_cluster <- (read_csv("path/output/t_clusters_LINKE.csv") %>% pull(label))[[1]]
CORPUS %>% kwic(paste0("[pcluster = '", gsub("\\|", "\\\\|", test_cluster),"']"),
                          s_attributes = c("party"),
                          cqp = TRUE)


keywords <- function(s_att, meta) {
  
  if (!(meta %in% names(KEYNESS_COMPARISONS))) {
    reference <- CORPUS
    included <- TRUE
  } else {
    reference <- CORPUS %>% partition(setNames(KEYNESS_COMPARISONS[[meta]], META_S_ATT[[meta]])) %>% enrich(p_attribute = "word")
    included <- FALSE
  }
  
  keywords <- CORPUS %>% partition(setNames(meta, s_att)) %>%
    enrich(p_attribute = "word") %>%
    features(reference, method = METHOD_KEYWORDS, included = included)
  
  keywords %>% saveRDS(
    file = paste0(
      "path/output/keywords_",
      meta,
      ".rds"
    )
  )
  keywords
}

collocations <- function(s_att, meta) {
  
  p_att <- META_P_ATT[[meta]]
  clusters <-
    CORPUS %>% partition(setNames(meta, s_att)) %>% enrich(p_attribute = p_att)
  
  collocs <- Cooccurrences(
    .Object = clusters,
    p_attribute = p_att,
    left = 5L,
    right = 5L,
  )
  if (METHOD_COLLOCS == "ll") {
    ll(collocs)
  } else if (METHOD_COLLOCS == "pmi") {
      pmi(collocs)
  } else if (METHOD_COLLOCS == "chisquare") {
    chisquare(collocs)
  } else if (METHOD_COLLOCS == "t_test") {
    t_test(collocs)
  } else {
    print("Please choose an available cooccurrence statistic!")
  }
  decode(collocs)
  collocs %>% saveRDS(
    file = paste0(
      "path/output/collocs_",
      meta,
      ".rds"
    )
  )
  collocs
}

cluster_frequencies <- function(s_att, meta) {
  
    cluster_freq <- hash()
    part <- CORPUS %>% partition(setNames(meta, s_att))
    for (cluster in names(META_CLUSTER_WORDS[[meta]])) {
      freq_df <-
        part %>% count(query = META_CLUSTER_WORDS[[meta]][[cluster]]) %>% tibble
      freq_df$freq <- round(freq_df$freq * 1000000, 2)
      cluster_freq[[cluster]] <- freq_df
    }
    cluster_freq %>% saveRDS(
      file = paste0(
        "path/output/cluster_freq_",
        meta,
        ".rds"
      )
    )
    cluster_freq
  }

fdgraph <- function(meta) {
  
  p_att <- META_P_ATT[["meta"]]
  clusters <- read_csv(
    paste0(
      "path/output/t_clusters_",
      meta,
      ".csv"
    )
  )
  column_a <- paste0("a_", META_P_ATT[[meta]])
  column_b <- paste0("b_", META_P_ATT[[meta]])
  collocs <- readRDS(
    file = paste0(
      "path/output/collocs_",
      meta,
      ".rds"
    )
  )@stat %>% tibble
  keywords <- readRDS(
    file = paste0(
      "path/output/keywords_",
      meta,
      ".rds"
    )
  )@stat
  cluster_freq <- readRDS(
    file = paste0(
      "path/output/cluster_freq_",
      meta,
      ".rds"
    )
  )
  t_collocs <- collocs %>%
    filter(column_a != "-" & column_b != "-") %>%
    mutate(
      a_id = map_chr(
        collocs[[column_a]],
        ~ clusters %>%
          filter(label == .x) %>%
          pull(cluster) %>%
          first() %>%
          as.character()
      ),
      b_id = map_chr(
        collocs[[column_b]],
        ~ clusters %>%
          filter(label == .x) %>%
          pull(cluster) %>%
          first() %>%
          as.character()
      )
    )
  nodes <-
    t_collocs %>%
    distinct(a_id, .keep_all = TRUE) %>%
    select(id = a_id,label = column_a) %>%
    filter(label != "-") %>%
    # create tooltip
    mutate(title = id %>%
             map(function(x) paste0(
                "<i>n = ",
                META_CLUSTER_WORDS[[meta]][[x]] %>% length %>% as.character,
                "</i><br><b>Keyness (LL) ",
                "n_key = ",
                keywords %>% filter(word %in% (META_CLUSTER_WORDS[[meta]][[x]]) &
                                                        ll >= 3.84) %>% nrow %>% as.character,
                " = ",
                (((keywords %>% filter(
                  word %in% (META_CLUSTER_WORDS[[meta]][[x]]) &
                    ll >= 3.84
                ) %>% nrow) / META_CLUSTER_WORDS[[meta]][[x]] %>% length
                ) * 100) %>% round(2) %>% as.character,
                " %</b><br>",
                keywords %>% filter(word %in% (META_CLUSTER_WORDS[[meta]][[x]]) &
                                                        ll >= 3.84) %>% head(15) %>% {
                                                          paste(.$word, round(.$ll, 2), sep = ": ")
                                                        } %>% toString() %>% str_replace_all(", ", "<br>"),
                "<br><br><b>Frequency</b><br>",
                cluster_freq[[x]] %>% arrange(desc(freq)) %>% select(query, freq) %>% head(15) %>% {
                  paste(.$query, .$freq, sep = ": ")
                } %>% toString() %>% str_replace_all(", ", "<br>")
                                  )
                )
          )
  
  edges <-
    t_collocs %>% select(from = a_id,
                         to = b_id,
                         value = ll)
  
  
  list(nodes, edges) %>% saveRDS(
    file = paste0(
      "path/output/fdgraph_",
      meta,
      ".rds"
    )
  )
  list(nodes, edges)
}

create_graph <- function(nodes, edges) {
  visNetwork(nodes, edges, height = "1250px", width = "100%") %>%
    visIgraphLayout(layout = "layout_with_mds") %>%
    visOptions(
      highlightNearest = list(
        enabled = T,
        degree = 1,
        hover = F,
        labelOnly = F
      ),
      nodesIdSelection = TRUE
    ) %>%
    visEdges(
      scaling = list(max = 5),
      color = list(opacity = 0.6),
      smooth = TRUE
    ) %>%
    visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(springConstant = 0.785),
      stabilization = F
    ) %>%
    visNodes(font = list(size = 14, strokeWidth = 4),
             size = 11) %>%
    visConfigure(enabled = T) %>%
    visInteraction(navigationButtons = TRUE, multiselect = TRUE) %>%
    visLayout(randomSeed = 1)
}

visualize_frame <- function(meta) {
  nodes_edges <- readRDS(
    file = paste0(
      "path/output/fdgraph_",
      meta,
      ".rds"
    )
  )
  
  
  t_clusters <- read_csv(
    paste(
      "path/output/t_clusters_",
      meta,
      ".csv",
      sep = ""
    )
  ) %>% filter(thematic_similarity >= T_SIM)
  
  nodes <-
    nodes_edges[[1]] %>% left_join(t_clusters, by = "label") %>%  filter(label %in% t_clusters$label) %>% mutate(id = id %>% as.integer) %>% arrange(id) %>% group_by(id) %>% mutate(id2 = cur_group_id()) %>% ungroup
  
  
  edges <-
    nodes_edges[[2]]  %>% mutate(
      from = from %>% as.integer,
      to = to %>% as.integer,
      id = from
    ) %>% left_join(nodes) %>% select(
      from = id2,
      to = to,
      value = value,
      id = to
    ) %>% left_join(nodes) %>% select(from = from,
                                      to = id2,
                                      weight = value) %>% arrange(from, to) %>% filter(from %in% nodes$id2 &
                                                                                         to %in% nodes$id2) %>% group_by(from) %>% slice_max(order_by = weight, n = 20) %>% ungroup
  
  nodes <-
    nodes %>% select(id = id2, label, title, n_words) %>% mutate(shape = case_when(n_words >= 200 ~ "star", n_words >= 100 ~ "diamond", TRUE ~ "dot"))
  
  
  g <- tbl_graph(nodes = nodes,
                 edges = edges,
                 directed = TRUE)
  
  g <- g %>% activate(nodes) %>% mutate(group = group_spinglass())
  
  nodes <- g %>% activate(nodes) %>% as_tibble
  edges <-
    g %>% activate(edges) %>% as_tibble %>% rename(value = weight)
  
  print(paste(nodes %>% pull(group) %>% unique %>% length, "groups"))
  
  list(nodes, edges) %>% saveRDS(
    paste0(
      "path/output/t_clusters_",
      meta,
      "_",
      gsub(".", "", T_SIM %>% as.character, fixed = T),
      "_visualisation.rds"
    )
  )
  
  
  create_graph(nodes, edges) %>%
    visSave(
      paste(
        "path/output/FDgraphs/FDgraph_",
        meta,
        "_",
        gsub(".", "", T_SIM %>% as.character, fixed = T),
        ".html",
        sep = ""
      )
    )
}

meta_keywords <- META_S_ATT %>% imap(keywords)
meta_collocs <- META_S_ATT %>% imap(collocations)
meta_cluster_word_freq <- META_S_ATT %>% imap(cluster_frequencies)
meta_fdgraph <- META_S_ATT %>% names %>% map(fdgraph)
META_S_ATT %>% names %>% map(visualize_frame)

### if you recognize recurring network communities in the generated FDgraphs, you can specify a label and a color for them
# you will have to cover all detected communities and specify one cluster of the community as an anchor
# this example will not work on your computer as you do not have identical clusters

community_cluster <- list()
community_cluster[["1990"]] <- list()
community_cluster[["1990"]][["1_Example1"]] <- "90_Chance|Gelegenheit"
community_cluster[["1990"]][["2_Example2"]] <- "90_Gewinn|Vorteil|Nutzen"
community_cluster[["1990"]][["3_Example3"]] <- "90_Familien|Jugendlichen|Eltern"
# community_cluster[["1990"]][["4_add_more_if_needed"]] <- ""

community_cluster[["2000"]] <- list()
community_cluster[["2000"]][["1_Example1"]] <- "00_Möglichkeit|Chance|Gelegenheit"
community_cluster[["2000"]][["2_Example2"]] <- "00_Vorteil|Preis|Effekt"
community_cluster[["2000"]][["3_Example3"]] <- "00_Eltern|Familien|Kindern"
# community_cluster[["2000"]][["4_add_more_if_needed"]] <- ""

community_cluster[["2010"]] <- list()
community_cluster[["2010"]][["1_Example1"]] <- "10_Möglichkeit|Chance|Gelegenheit"
community_cluster[["2010"]][["2_Example2"]] <- "10_Bonus|Abschlag|Rabatt"
community_cluster[["2010"]][["3_Example3"]] <- "10_Kinder|Kindern"
# community_cluster[["2010"]][["4_add_more_if_needed"]] <- ""

fix_groups <-  function(meta){
  
  nodes_edges <- readRDS(
    paste0(
      "path/output/t_clusters_",
      meta,
      "_",
      gsub(".", "", T_SIM %>% as.character, fixed = T),
      "_visualisation.rds"
    )
  )
  nodes <- nodes_edges[[1]]
  edges <- nodes_edges[[2]]

  example1 <-
    ifelse(group_cluster[[meta]][["1_Example1"]] != "",
           (nodes %>% filter(label == group_cluster[[meta]][["1_Example1"]]))[["group"]] %>% as.character,
           "-1")
  
  example2 <-
    ifelse(group_cluster[[meta]][["2_Example2"]] != "",
           (nodes %>% filter(label == group_cluster[[meta]][["2_Example2"]]))[["group"]] %>% as.character,
           "-1")
  
  example3 <-
    ifelse(group_cluster[[meta]][["3_Example3"]] != "",
           (nodes %>% filter(label == group_cluster[[meta]][["3_Example3"]]))[["group"]] %>% as.character,
           "-1")
  
  #example4 <- ... add more if needed
  
  
  tryCatch({
    nodes <-
      nodes %>% mutate(group = if_else(group == example1, "1_Example1", group %>% as.character))
    nodes <-
      nodes %>% mutate(group = if_else(group == example2, "2_Example2", group %>% as.character))
    nodes <-
      nodes %>% mutate(group = if_else(group == example3, "3_Example3", group %>% as.character))
    # nodes <- ... add more if needed
  }, error = function(e) {})
  
  nodes <- nodes %>% arrange(group, label)
  
  list(nodes, edges) %>% saveRDS(
    paste0(
      "path/output/t_clusters_",
      meta,
      "_",
      gsub(".", "", e_sim %>% as.character, fixed = T),
      "_visualisation_fixed_groups.rds"
    )
  )
  
  visNetwork(nodes, edges, height = "1250px", width = "100%") %>%
    visIgraphLayout("layout_with_mds") %>%
    visOptions(
      highlightNearest = list(
        enabled = T,
        degree = 1,
        hover = F,
        labelOnly = F
      ),
      nodesIdSelection = TRUE
    ) %>%
    visEdges(
      scaling = list(max = 5),
      color = list(opacity = 0.6),
      smooth = TRUE
    ) %>%
    visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(springConstant = 0.785),
      stabilization = F
    ) %>%
    visNodes(font = list(size = 14, strokeWidth = 4),
             size = 11) %>%
    visGroups(
      groupname = "1_Example1",
      color = list(
        border = "#2B7CE9",
        background = "#97C2FC",
        highlight = list(border = "#2B7CE9", background = "#D2E5FF"),
        hover = list(border = "#2B7CE9", background = "#D2E5FF")
      )
    ) %>%
    visGroups(
      groupname = "2_Example2",
      color = list(
        border = "#FFA500",
        background = "#FFFF00",
        highlight = list(border = "#FFA500", background = "#FFFFA3"),
        hover = list(border = "#FFA500", background = "#FFFFA3")
      )
    ) %>%
    visGroups(
      groupname = "3_Example3",
      color = list(
        border = "#FA0A10",
        background = "#FB7E81",
        highlight = list(border = "#FA0A10", background = "#FFAFB1"),
        hover = list(border = "#FA0A10", background = "#FFAFB1")
      )
    )  %>%
    visLegend(width = 0.08, position = "right") %>%
    visConfigure(enabled = T) %>%
    visInteraction(navigationButtons = TRUE, multiselect = TRUE) %>%
    visLayout(randomSeed = 1) %>%
    visSave(
      paste(
        "path/output/FDgraphs/FDgraph_",
        meta,
        "_",
        gsub(".", "", e_sim %>% as.character, fixed = T),
        "_fixed_groups.html",
        sep = ""
      )
    )
}
