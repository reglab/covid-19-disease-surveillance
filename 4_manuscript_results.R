library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(forcats)
library(lubridate)
library(jsonlite)
library(phangorn)
library(xtable)
set.seed(1)

cat("=== GETTING PROPOSED EDGES... ===", "\n")
source("3_get_proposed_edges.R")
from_labels <- all_proposed_edges$from
to_labels <- all_proposed_edges$to
num_pairs <- length(from_labels)

cat("=== READING GISAID & TERRA DATA ===", "\n")
strains <- read.phyDat("[REDACTED]", format="fasta") # cannot be shared
metadata <- read_csv("[REDACTED]") # cannot be shared
coll_dates <- metadata %>% 
  select(strain, collection_date)
cc_to_strains <- read_csv("[REDACTED]") # cannot be shared
strains <- strains[names(strains) %in% cc_to_strains$strain]
sources <- as.data.frame(names(strains)) %>% dplyr::rename(strain = 1) %>% 
  left_join(metadata %>% select(strain, source), by="strain")

cat("=== RESAMPLING CASES... ===", "\n")
cc_case_resamples <- lapply(1:100, function(x){sample(calconnect_dataset$CC_RECORD_NUM, length(calconnect_dataset$CC_RECORD_NUM), replace=T)})
cc_resample_strain_counts <- lapply(cc_case_resamples, function(cases) {
  cases %>% as_tibble(.name_repair="unique") %>% 
    dplyr::rename("CC_RECORD_NUM" = 1) %>% 
    left_join(cc_to_strains, by=c("CC_RECORD_NUM")) %>%
    filter(!is.na(strain)) %>% 
    group_by(strain) %>% count()
})

cat("=== GETTING BOOTSTRAPPED DISTANCES FROM TREES ===")
trees <- c()
for (file_name in list.files("[REDACTED]")) {
  new_trees <- read.tree(paste0("[REDACTED]", file_name)) # cannot be shared
  trees <- c(trees, new_trees)
}
cat("Read", length(trees), "trees.\n")
bootstrap_dists <- all_proposed_edges
i <- 0
for (tree in trees) {
  i <- i + 1
  cat("Tree", i, "\n")
  cophenetic_result <- cophenetic.phylo(tree)
  relevant_dists <- vapply(1:num_pairs, FUN=function(i) {
    cophenetic_result[from_labels[i], to_labels[i]]
  }, numeric(1))
  bootstrap_dists <- bootstrap_dists %>% add_column(relevant_dists, .name_repair="universal")
}
names(bootstrap_dists) <- c(names(all_proposed_edges), paste0("iter", 1:length(trees)))

cat("=== REWEIGHTING PROPOSED EDGES... ===", "\n")
weighted_proposed_edges <- lapply(cc_resample_strain_counts, function(strain_counts) {
  all_proposed_edges %>% 
    left_join(strain_counts, by=c("from" = "strain")) %>% 
    left_join(strain_counts, by=c("to" = "strain")) %>% 
    replace_na(list(n.x = 0, n.y = 0)) %>% 
    mutate(weight = n.x * n.y)
})

threshold = 1 / 29903
v4_performance_matrix_threshold1 <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
threshold1_means <- rowMeans(v4_performance_matrix_threshold1, na.rm=T)
threshold1_quantiles <- apply(v4_performance_matrix_threshold1, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance_threshold1 <- bind_cols(method=list_of_methods, sim=threshold1_means, ci_lower=threshold1_quantiles[, 1], ci_upper=threshold1_quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges$edgename)
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))

threshold = 2 / 29903
v4_performance_matrix <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
means <- rowMeans(v4_performance_matrix, na.rm=T)
quantiles <- apply(v4_performance_matrix, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance <- bind_cols(method=list_of_methods, sim=means, ci_lower=quantiles[, 1], ci_upper=quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges$edgename)
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))
saveRDS(v4_performance, "v4_performance.rds")

# robustness check, changing threshold----
threshold = 3 / 29903
v4_performance_matrix_threshold3 <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
threshold3_means <- rowMeans(v4_performance_matrix_threshold3, na.rm=T)
threshold3_quantiles <- apply(v4_performance_matrix_threshold3, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance_threshold3 <- bind_cols(method=list_of_methods, sim=threshold3_means, ci_lower=threshold3_quantiles[, 1], ci_upper=threshold3_quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges$edgename)
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))

threshold = 4 / 29903
v4_performance_matrix_threshold4 <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
threshold4_means <- rowMeans(v4_performance_matrix_threshold4, na.rm=T)
threshold4_quantiles <- apply(v4_performance_matrix_threshold4, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance_threshold4 <- bind_cols(method=list_of_methods, sim=threshold4_means, ci_lower=threshold4_quantiles[, 1], ci_upper=threshold4_quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges$edgename)
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))

threshold = 5 / 29903
v4_performance_matrix_threshold5 <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
threshold_means <- rowMeans(v4_performance_matrix_threshold5, na.rm=T)
threshold_quantiles <- apply(v4_performance_matrix_threshold5, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance_threshold5 <- bind_cols(method=list_of_methods, sim=threshold_means, ci_lower=threshold_quantiles[, 1], ci_upper=threshold_quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges$edgename)
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))

# Check how many of all bootstrapped distances are less than----
test <- 1:100 %>% map(function(x){
  
  iter <- bootstrap_dists %>%
    pull(paste0("iter",x))
  
  sum(iter < 3/29903)
  
  
}) %>% unlist() %>% sum()
test/nrow(bootstrap_dists)/100

# Alternative measures comparison----

cat("=== VERSION 1: SNP DISTANCES, BINOMIAL INTERVAL ===", "\n")
snp_dists <- read_tsv("[REDACTED]") %>% # cannot be shared
  pivot_longer(cols=2:11984, names_to="seq2", values_to="dist") %>%
  dplyr::rename(seq1 = 1) %>%
  mutate(edgename = paste0(seq1, "->", seq2)) %>%
  select(edgename, dist)

snp_ground_truth <- snp_dists %>%
  filter(dist < 2)

v1_performance <- sapply(list_of_methods, function(method) {
  
  # Sequenced edges proposed by method
  proposed_edges <- get(method)$all_sequenced_edges$edgename
  
  # Sequenced edges proposed by method, confirmed similar by tree
  confirmed_edges <- intersect(proposed_edges, snp_ground_truth$edgename)
  
  similarity <- length(confirmed_edges) / length(proposed_edges)
  ci <- binom.test(x = length(confirmed_edges), n = length(proposed_edges), conf.level=0.95)$conf.int
  return(c(similarity, ci[1], ci[2]))
}, simplify="matrix") %>% t() %>% as.data.frame() %>%
  dplyr::rename(sim=1, ci_lower=2, ci_upper=3) %>%
  mutate(width = abs(ci_lower - ci_upper)) %>%
  rownames_to_column(var="method")

cat("=== VERSION 2: SNP DISTANCES, CASE-BOOTSTRAP INTERVAL ===", "\n")
v2_performance_matrix <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>%
    mutate(edgename = paste0(from, "->", to)) %>%
    select(edgename, weight)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, snp_ground_truth$edgename)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
}, simplify="matrix")
means <- rowMeans(v2_performance_matrix, na.rm=T)
quantiles <- apply(v2_performance_matrix, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v2_performance <- bind_cols(method=list_of_methods, sim=means, ci_lower=quantiles[, 1], ci_upper=quantiles[, 3]) %>%
  mutate(width = abs(ci_upper - ci_lower))

cat("=== VERSION 3: BOOTSTRAP DISTANCES, PHYLO-UNCERTAINTY INTERVAL ===", "\n")
v3_performance_matrix <- sapply(1:100, function(i) {
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>%
    select(from, to, dists = one_of(paste0("iter", i))) %>%
    mutate(edgename = paste0(from, "->", to)) %>%
    filter(dists < threshold) %>%
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges$edgename
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Similarity quotient
    similarity <- length(confirmed_edges) / length(proposed_edges)
    return(similarity)
  })
  
}, simplify="matrix")
means <- rowMeans(v3_performance_matrix, na.rm=T)
quantiles <- apply(v3_performance_matrix, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v3_performance <- bind_cols(method=list_of_methods, sim=means, ci_lower=quantiles[, 1], ci_upper=quantiles[, 3]) %>%
  mutate(width = abs(ci_upper - ci_lower))

# Set ordering for methods for downstream analysis
v4_performance <- readRDS("v4_performance.rds")

levels <- v4_performance %>% arrange(sim) %>% pull(method)

methods_to_show <- v4_performance %>% arrange(desc(sim)) %>% filter(method %in% c(
  "same_address", "location_history",
  "work_exposure", "same_employer_name", "same_supervisor", 
  "occ_location", "school_exposure", "same_messy_school",
  "all_congregate_exposure", "case_contact_edges", "satscan_precise",
  "ct_excluding_same_address", "satscan_excluding_same_address", 
  "ltcf_exposure", "jails_exposure")) %>% pull(method)

map_methods <- function(methods) {
  vapply(methods, function(x) {
    new_name <- switch(x,
                       "same_cbg" = "Same CBG",
                       "same_address" = "Same Address",
                       "parcel" = "Same Parcel",
                       "location_history" = "Location History",
                       "work_exposure" = "Workplace Reporting",
                       "same_messy_school" = "Same School",
                       "school_exposure" = "School Reporting",
                       "same_supervisor" = "Same Supervisor",
                       "same_employer_name" = "Same Employer Name",
                       "occ_location" = "Same Occupation Location",
                       "case_contact_edges" = "Contact Tracing",
                       "satscan_precise" = "SaTScan",
                       "satscan_excluding_same_address" = "SaTScan, excl. Same Address",
                       "ct_excluding_same_address" = "CT, excl. Same Address",
                       "ltcf_exposure" = "LTCF Reporting",
                       "all_congregate_exposure" = "All Mandated Reporting",
                       "jails_exposure" = "Jail Reporting",
                       "unhoused_exposure" = "Unhoused Reporting",
                       "UNKNOWN"
    )
    if (new_name == "UNKNOWN") {return(x)}
    return(new_name)
  }, character(1)) %>% unname()
}

### GRAPHS ###

# MAIN FIGURE: Comparing surveillance methodologies ----
# Visual Part
pdf("[REDACTED]", width=5.5, height=8)
v4_performance %>% 
  filter(method %in% methods_to_show) %>% 
  arrange(sim) %>% 
  ggplot(aes(y = factor(method, levels=levels), x = sim, xmin = ci_lower, xmax = ci_upper)) + 
  geom_point(aes(size=sample_size), color="#D3D3D3") + 
  geom_errorbar(width=0.2, color="#4E79A7") +
  geom_point(size=1, color="#4E79A7") + 
  geom_hline(yintercept = 1:(length(methods_to_show) - 1) + 0.5, color="gray") +
  labs(
    title="",
    x = "",
    y = "",
    size = "Sequenced\nProposed Links"
  ) +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.line.y = element_blank(),
  ) + 
  scale_size_continuous(breaks=c(100, 500, 2000, 8000, 12000), range = c(5,18), guide="none")
dev.off()

v4_performance %>% 
  filter(method %in% methods_to_show) %>% 
  arrange(desc(sim)) %>% 
  transmute(
    `Surveillance Method` = map_methods(method),
    `Informational Value` = format(sim, digits=2),
    `95% CI` = paste0("(", format(ci_lower, digits=1), ", ", format(ci_upper, digits=2), ")"),
    `# Sequenced Links` = sample_size
  ) %>% xtable() %>% print.xtable(type="html", include.rownames=F, file="[REDACTED]")


# Fig S1: Full Results & Overlap ----
overlap <- function(t1, t2, edgetype) {
  if (edgetype == "unsequenced") {
    edges1 <- get(t1)$all_unsequenced_edges$edgename
    edges2 <- get(t2)$all_unsequenced_edges$edgename
  } else if (edgetype == "sequenced") {
    edges1 <- get(t1)$all_possible_edges$edgename
    edges2 <- get(t2)$all_possible_edges$edgename
  } else if (edgetype == "ground truth") {
    edges1 <- get(t1)$all_possible_edges$edgename
    edges1 <- edges1[edges1 %in% genomic_edges$edgename]
    edges2 <- get(t2)$all_possible_edges$edgename
    edges2 <- edges2[edges2 %in% genomic_edges$edgename]
  }
  
  jaccard <- length(base::intersect(edges1, edges2)) / length(base::union(edges1, edges2))
  cp <- length(unique(base::intersect(edges1, edges2))) / length(unique(edges2))
  tribble(~method1, ~method2, ~jaccard, ~p_1_given_2,
          t1, t2, jaccard, cp)
}

calculate_overlap_matrix <- function(methods, edgetype) {
  res <- data.frame()
  for (m in methods) {
    for (n in methods) {
      res <- bind_rows(res, overlap(m, n, edgetype))
    }
  }
  res
}

overlap_df <- calculate_overlap_matrix(methods_to_show, "unsequenced")

pdf("[REDACTED]", width=14, height=10)
overlap_df %>% 
  mutate(
    method1 = map_methods(method1),
    method2 = map_methods(method2),
    p_1_given_2 = if_else(method1 == method2, as.numeric(NA), p_1_given_2)
    
  ) %>% 
  ggplot() +
  geom_tile(aes(y = fct_rev(factor(method1, levels=map_methods(methods_to_show))), 
                x = factor(method2, levels=map_methods(methods_to_show)),
                fill = p_1_given_2)) +
  geom_text(aes(y = fct_rev(factor(method1, levels=map_methods(methods_to_show))), 
                x = factor(method2, levels=map_methods(methods_to_show)),
                label = round(p_1_given_2, 2))) +
  scale_fill_distiller(type="seq", palette="YlOrRd", direction=1) +
  guides(x = guide_axis(angle = 60)) +
  scale_x_discrete(position = "top", labels=map_methods(methods_to_show)) +
  scale_y_discrete(labels = rev(map_methods(methods_to_show))) +
  labs(x = "Strategy 2", y = "Strategy 1", 
       fill = "P(S1 | S2)") +
  theme_minimal(base_size = 20)
dev.off()

# Fig S2:  P-values ----
p_vals <- sapply(methods_to_show, function(method1) {
  sapply(methods_to_show, function(method2) {
    mean(v4_performance_matrix[method1,] > v4_performance_matrix[method2,])
  })
}) %>% as.data.frame() %>% 
  rownames_to_column(var="method") %>%
  pivot_longer(cols=-method, names_to = "method2", values_to="p")

jpeg("[REDACTED]", width=14, height=10)
p_vals %>% 
  rowwise() %>% 
  mutate(
    method1 = map_methods(method),
    method2 = map_methods(method2),
    p = case_when(
      method1 == method2 ~ as.numeric(NA),
      T ~ p
    ),
    label = case_when(
      !is.na(p) ~ sprintf("%.2f", p),
      T ~ ""
    )
  ) %>% 
  ungroup() %>% 
  ggplot() +
  geom_tile(aes(x = factor(method1, levels=map_methods(methods_to_show)), 
                y = fct_rev(factor(method2, levels=map_methods(methods_to_show))),
                fill = p)) +
  geom_text(aes(factor(method1, levels = map_methods(methods_to_show)), 
                y = fct_rev(factor(method2, levels = map_methods(methods_to_show))),
                label = label)) +
  scale_fill_distiller(type="div", palette="RdBu", na.value="grey50", direction=1) +
  guides(x = guide_axis(angle = 60)) +
  scale_x_discrete(position = "top", labels = map_methods(methods_to_show)) +
  scale_y_discrete(labels = rev(map_methods(methods_to_show))) +
  labs(x = "Strategy 2", y = "Strategy 1",  
       fill = "P(S1 > S2)") +
  theme_minimal(base_size = 20)
dev.off()

# APPENDIX I: GISAID vs. Terra Robustness Check
groups = split(sources, sources$source)
terra_strains <- groups$`Terra/PHAGE` %>% pull(strain)
gisaid_strains <- groups$`GISAID/Fulgent` %>% pull(strain)

### Terra only
v4_performance_matrix_terra <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges %>% 
      filter(from %in% terra_strains) %>% pull(edgename)
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
terra_means <- rowMeans(v4_performance_matrix_terra, na.rm=T)
terra_quantiles <- apply(v4_performance_matrix_terra, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance_terra <- bind_cols(
  method=list_of_methods, 
  sim=terra_means, 
  ci_lower=terra_quantiles[, 1], 
  ci_upper=terra_quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges %>% 
                           filter(from %in% terra_strains) %>% pull(edgename))
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))

### GISAID only
v4_performance_matrix_gisaid <- sapply(1:100, function(i) {
  # Based on case resampling
  edge_weights_i <- weighted_proposed_edges[[i]] %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    select(edgename, weight)
  
  # Based on bootstrap resampling
  ground_truth_i <- bootstrap_dists %>% 
    select(from, to, dists = one_of(paste0("iter", i))) %>% 
    mutate(edgename = paste0(from, "->", to)) %>% 
    filter(dists < threshold) %>% 
    pull(edgename)
  
  sapply(list_of_methods, function(method) {
    # Sequenced edges proposed by method
    proposed_edges <- get(method)$all_sequenced_edges %>% 
      filter(from %in% gisaid_strains) %>% pull(edgename)
    
    # Sequenced edges proposed by method, confirmed similar by tree
    confirmed_edges <- intersect(proposed_edges, ground_truth_i)
    
    # Weighted count of  proposed edges
    proposed_sum <- edge_weights_i %>%
      filter(edgename %in% proposed_edges) %>%
      pull(weight) %>% sum()
    
    # Weighted count of confirmed edges
    confirmed_sum <- edge_weights_i %>%
      filter(edgename %in% confirmed_edges) %>%
      pull(weight) %>% sum()
    
    # Similarity quotient
    similarity <- confirmed_sum / proposed_sum
    return(similarity)
  })
  
}, simplify="matrix")
gisaid_means <- rowMeans(v4_performance_matrix_gisaid, na.rm=T)
gisaid_quantiles <- apply(v4_performance_matrix_gisaid, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T) %>% t()
v4_performance_gisaid <- bind_cols(
  method=list_of_methods, 
  sim=gisaid_means, 
  ci_lower=gisaid_quantiles[, 1], 
  ci_upper=gisaid_quantiles[, 3]) %>% 
  rowwise() %>% 
  mutate(
    width = abs(ci_upper - ci_lower),
    sample_size = length(get(method)$all_sequenced_edges %>% 
                           filter(from %in% gisaid_strains) %>% pull(edgename))
  ) %>%
  ungroup() %>% 
  arrange(desc(sim))

# Fig S3: Comparing performance by method (as dot plot) ----
combined_gisaid_terra <- bind_rows(
  v4_performance_terra %>% mutate(source="Terra"),
  v4_performance_gisaid %>% mutate(source="GISAID")
) %>% 
  filter(
    !(method %in% c("all_congregate_exposure", "jails_exposure", "ltcf_exposure")),
    method %in% methods_to_show
  )

pdf("[REDACTED]", width=12, height=6)
combined_gisaid_terra %>% 
  filter(method %in% methods_to_show) %>% 
  ggplot(aes(x = sim, y = factor(map_methods(method), levels=(map_methods(levels))), 
             xmax = ci_upper, xmin=ci_lower, color=factor(source))) + 
  geom_point(size=2.5, position=position_dodge2(width=0.9)) +
  geom_hline(yintercept = 1:length(methods_to_show) + 0.5, color="gray") +
  scale_color_discrete(guide=guide_legend(reverse=T)) +
  theme_minimal() +
  labs(
    y = "Surveillance Strategy", x = "Informational Value", color="Genome Source")
dev.off()

# Fig S4: Correlation of performance of two methods ----
corr_df <- v4_performance_terra %>% 
  filter(
    !(method %in% c("all_congregate_exposure", "jails_exposure", "ltcf_exposure")),
    method %in% methods_to_show
  ) %>% select(method, sim) %>% 
  left_join(v4_performance_gisaid %>% 
              filter(
                !(method %in% c("all_congregate_exposure", "jails_exposure", "ltcf_exposure")),
                method %in% methods_to_show
              ) %>% select(method, sim), by="method") %>% 
  filter(method %in% methods_to_show)

pdf("[REDACTED]", width=12, height=6)
corr_df %>% 
  ggplot(aes(x = sim.x, y = sim.y)) +
  geom_point(color="#4E79A7") +
  geom_smooth(method='lm', color="#F28E2B", fill="gray90", formula=y ~ x) +
  theme_minimal() +
  labs(
    y = "GISAID Informational Value", x = "Terra Informational Value")
dev.off()

cat("Pearson Correlation: ", cor(corr_df$sim.x, corr_df$sim.y, method="pearson"), "\n")
cat("Spearman Correlation: ", cor(corr_df$sim.x, corr_df$sim.y, method="spearman"), "\n")

# Fig S5: SaTScan Hyperparameter Tuning Chart ----
results_dir <- "[REDACTED]"
json_files <- list.files(results_dir) %>% 
  as.data.frame() %>% dplyr::rename("file_name" = ".") %>% 
  filter(str_detect(file_name, "json")) %>% 
  separate(file_name, into = c("x", "time", "y", "radius", "z", "q", "idx"), sep="-", remove=F) %>%
  transmute(file_name = file_name, time = as.numeric(time), 
            radius = as.numeric(str_replace(radius, "_", ".")),
            idx = as.numeric(idx %>% str_replace("\\.json", ""))) %>% 
  arrange(desc(idx))
cc_with_strain <- calconnect_dataset %>% filter(!is.na(strain))
with_household_and_strain <- calconnect_dataset %>% filter(!is.na(strain), !is.na(HOUSEHOLD))
household_edges <- with_household_and_strain %>% 
  left_join(with_household_and_strain, by = c("HOUSEHOLD")) %>% 
  filter(strain.x < strain.y) %>% 
  select(from = strain.x, to = strain.y,
  ) %>% mutate(edgename = paste0(from, "->", to))

with_melissa_and_strain <- calconnect_dataset %>% 
  filter(!is.na(strain), !is.na(MELISSA_ID), !is.na(PRECISION)) %>% 
  filter(PRECISION != "NEARBY ADDRESS")
same_address_edges <- with_melissa_and_strain %>% 
  left_join(with_melissa_and_strain, by = c("MELISSA_ID")) %>% 
  filter(strain.x < strain.y) %>% 
  select(from = strain.x, to = strain.y) %>% 
  mutate(edgename = paste0(from, "->", to))

get_satscan_edges <- function(results_dir, json_file, calconnect) {
  # Read in cluster information
  data_path <- paste0(results_dir, json_file)
  clusters <- read_json(path=data_path, simplifyDataFrame=T) %>% 
    as_tibble() %>% 
    unnest_longer(col=points) %>% 
    unnest_wider(col=points) %>%
    dplyr::rename(`point_lng` = `...1`, `point_lat` = `...2`)
  
  # Join to CalCONNECT to recover individual records for each cluster
  with_strain <- clusters %>% 
    mutate(
      start_date = as.Date(start_date, format="%Y/%m/%d"),
      end_date = as.Date(end_date, format="%Y/%m/%d")) %>% 
    left_join(calconnect, by=c("point_lat" = "LAT", "point_lng" = "LON")) %>% 
    filter(EPISODE_DATE >= start_date, EPISODE_DATE <= end_date) %>%
    filter(!is.na(strain))
  
  all_sequenced_edges <- with_strain %>% 
    left_join(with_strain, by=c("id")) %>% 
    filter(strain.x < strain.y) %>% 
    select(from = strain.x, to = strain.y,
           from_episode_date = EPISODE_DATE.x,
           to_episode_date = EPISODE_DATE.y
    ) %>% 
    mutate(
      edgename = paste0(from, "->", to),
      time_diff = as.numeric(abs(from_episode_date - to_episode_date))
    ) %>% filter(time_diff < 15) %>% 
    group_by(edgename) %>% 
    dplyr::slice(1) %>% ungroup()
  
  # Calculate possible/actual edges
  proposed_edges <- all_sequenced_edges$edgename
  ground_truth <- intersect(snp_ground_truth$edgename, proposed_edges)
  
  inter_household_proposed <- setdiff(proposed_edges, household_edges$edgename)
  inter_household_ground_truth <- intersect(snp_ground_truth$edgename, inter_household_proposed)
  
  inter_address_proposed <- setdiff(proposed_edges, same_address_edges$edgename)
  inter_address_ground_truth <- intersect(snp_ground_truth$edgename, inter_address_proposed)
  
  list(json_file, proposed_edges, ground_truth, 
       inter_household_proposed, inter_household_ground_truth,
       inter_address_proposed, inter_address_ground_truth)
}

tuning_results <- mclapply(json_files$file_name, function(file_name) {
  get_satscan_edges(results_dir, file_name, cc_with_strain)
}, mc.cores=10) 

tuning_df <- sapply(tuning_results, function(item) {
  c(item[[1]], length(unique(item[[2]])), length(unique(item[[3]])), 
    length(unique(item[[4]])), length(unique(item[[5]])),
    length(unique(item[[6]])), length(unique(item[[7]])))
}) %>% t() %>% as.data.frame() %>% type.convert() %>% as_tibble()
names(tuning_df) <- c(
  "file_name", "proposed_edges", "ground_truth", 
  "inter_household_proposed", "inter_household_ground_truth",
  "inter_address_proposed", "inter_address_ground_truth")

wider <- json_files %>% 
  left_join(tuning_df, by=c("file_name"))

with_similarity <- wider %>% 
  mutate(
    all_incorrect_edges = proposed_edges - ground_truth,
    inter_household_incorrect = inter_household_proposed - inter_household_ground_truth,
    inter_address_incorrect = inter_address_proposed - inter_address_ground_truth,
    radius_quintile = ntile(radius, 5)
  )

pdf("[REDACTED]", width=13, height=8)
ggplot(with_similarity) +
  geom_point(aes(x = inter_address_ground_truth, y = inter_address_incorrect)) +
  geom_smooth(aes(x = inter_address_ground_truth, y = inter_address_incorrect), method="lm") +
  facet_grid(cols=vars(radius_quintile), rows=vars(time),
             labeller=labeller(radius_quintile = function(x) {paste("Radius Quintile", x)},
                               time = function(x) {paste(x, "days")})) +
  labs(x = "Inter-Address True Positive Edges", y = "Inter-Address False Positive Edges")
dev.off()
write_csv(with_similarity, "randomized-satscan-hyperparam-results.csv")

# Fig S6: SaTScan Calibration Results Table ----

# Visual Part
satscan_methods <- c("case_contact_edges", "satscan_precise", "satscan_loose")
pdf("[REDACTED]", width=5.5, height=2.5)
v4_performance %>% 
  filter(method %in% satscan_methods) %>% 
  arrange(sim) %>% 
  ggplot(aes(y = factor(method, levels=rev(satscan_methods)), x = sim, xmin = ci_lower, xmax = ci_upper)) + 
  geom_point(aes(size=sample_size), color="#D3D3D3") + 
  geom_errorbar(width=0.2, color="#4E79A7") +
  geom_point(size=1, color="#4E79A7") + 
  geom_hline(yintercept = 1:(length(satscan_methods) - 1) + 0.5, color="gray") +
  labs(title="", x = "", y = "", size = "Sequenced\nProposed Links") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.line.y = element_blank()) + 
  scale_size_continuous(breaks=c(100, 500, 2000, 8000, 12000), range = c(5,20), guide="none")
dev.off()

# HTML Table Part
v4_performance %>% 
  filter(method %in% satscan_methods) %>% 
  arrange(desc(sim)) %>% 
  transmute(
    `Surveillance Method` = map_methods(method),
    `Informational Value` = format(sim, digits=2),
    `95% CI` = paste0("(", format(ci_lower, digits=1), ", ", format(ci_upper, digits=2), ")"),
    `# Sequenced Links` = sample_size
  ) %>% xtable() %>% print.xtable(type="html", include.rownames=F, file="[REDACTED]")

# Fig S7: Case Count and Sequence Availability ----
by_week <- calconnect_dataset %>% 
  group_by(EPISODE_WEEK) %>% 
  summarize(n = n(), n_with_strain = sum(!is.na(strain)), seq_rate = n_with_strain / n) %>% 
  pivot_longer(cols = c(seq_rate, n_with_strain, n), names_to = "type", values_to = "value") %>% 
  mutate(
    type = case_when(
      type == "seq_rate" ~ "Sequencing Rate", type == "n_with_strain" ~ "Count With High-\nQuality Sequences",
      T ~ "Count of Cases")
  )

pdf("[REDACTED]", width=17, height=12)
ggplot(by_week) + 
  geom_col(aes(EPISODE_WEEK, value, fill=type)) +
  scale_x_date(date_breaks = "months" , date_labels = "%b '%y") +
  guides(x = guide_axis(angle = 45)) +
  facet_grid(type ~ ., scales="free_y", switch="both") +
  scale_fill_tableau(guide="none") +
  labs(title = "", x = "COVID-19 Episode Date", y = "") +
  theme_minimal(base_size = 16) + theme(strip.placement = "outside")
dev.off()

# Balance Tables
n_all_cases <- nrow(calconnect_dataset)
n_all_sequences <- sum(!is.na(calconnect_dataset$strain))
all_seq_rate <- n_all_sequences / n_all_cases
all_ci <- binom.test(x = n_all_sequences, n = n_all_cases, conf.level = 0.95)$conf.int

# Table S1: Balance Table by Demographic Group ---- 
demog_table <- tribble(~group, ~n_cases, ~n_sequences, ~seq_rate, ~ci_lower, ~ci_upper,
                       "All", n_all_cases, n_all_sequences, all_seq_rate, all_ci[[1]], all_ci[[2]])

# By Gender
calconnect_dataset <- calconnect_dataset %>% 
  mutate(
    GENDER_SIMPLE = case_when(
      str_detect(GENDER, "Female") ~ "Female",
      str_detect(GENDER, "Male") | str_detect(GENDER, "male") ~ "Male",
      T ~ "Other or Unknown"
    ))
for (gender in c("Female", "Male", "Other or Unknown")) {
  subset <- calconnect_dataset %>% filter(GENDER_SIMPLE == gender)
  cases <- nrow(subset)
  sequences <- sum(!is.na(subset$strain))
  rate = sequences / cases
  ci <- binom.test(x = sequences, n = cases, conf.level = 0.95)$conf.int
  trow <- tribble(~group, ~n_cases, ~n_sequences, ~seq_rate, ~ci_lower, ~ci_upper,
                  gender, cases, sequences, rate, ci[[1]], ci[[2]])
  demog_table <- bind_rows(demog_table, trow)
}

# By Race
calconnect_dataset <- calconnect_dataset %>% 
  mutate(
    RACE_SIMPLE = case_when(
      RACE == "White" ~ "White",
      RACE == "Black or African American" ~ "Black",
      RACE == "Asian" ~ "Asian",
      RACE == "Unknown" ~ "Unknown Race",
      str_detect(RACE, "Hawaiian") ~ "Hawaiian / Pacific Islander",
      str_detect(RACE, "Indian") ~ "American Indian / Alaska Native",
      T ~ "Other / Two or More Races"
    ))

for (race in unique(calconnect_dataset$RACE_SIMPLE)) {
  subset <- calconnect_dataset %>% filter(RACE_SIMPLE == race)
  cases <- nrow(subset)
  sequences <- sum(!is.na(subset$strain))
  rate = sequences / cases
  ci <- binom.test(x = sequences, n = cases, conf.level = 0.95)$conf.int
  trow <- tribble(~group, ~n_cases, ~n_sequences, ~seq_rate, ~ci_lower, ~ci_upper,
                  race, cases, sequences, rate, ci[[1]], ci[[2]])
  demog_table <- bind_rows(demog_table, trow)
}

# By Ethnicity
calconnect_dataset <- calconnect_dataset %>% 
  mutate(
    ETHNICITY_SIMPLE = case_when(
      ETHNICITY == "Hispanic or Latino" ~ "Hispanic/Latino",
      ETHNICITY == "Not Hispanic or Latino" ~ "Not Hispanic/Latino",
      T ~ "Other/Unknown Ethnicity"
    ))

for (ethn in unique(calconnect_dataset$ETHNICITY_SIMPLE)) {
  subset <- calconnect_dataset %>% filter(ETHNICITY_SIMPLE == ethn)
  cases <- nrow(subset)
  sequences <- sum(!is.na(subset$strain))
  rate = sequences / cases
  ci <- binom.test(x = sequences, n = cases, conf.level = 0.95)$conf.int
  trow <- tribble(~group, ~n_cases, ~n_sequences, ~seq_rate, ~ci_lower, ~ci_upper,
                  ethn, cases, sequences, rate, ci[[1]], ci[[2]])
  demog_table <- bind_rows(demog_table, trow)
}
write_csv(demog_table, "[REDACTED]")

# Table S2: Balance Table by Surveillance Method ----
all_pct_hispanic = calconnect_dataset %>% filter(ETHNICITY_SIMPLE == "Hispanic/Latino") %>% nrow() / n_all_cases
seq_rate_by_method <- tribble(~method, ~num_cases, ~num_sequences, ~seq_rate, ~ci_lower, ~ci_upper, ~pct_hisp,
                              "all", n_all_cases, n_all_sequences, all_seq_rate, all_ci[[1]], all_ci[[2]], all_pct_hispanic)
for (method in levels[levels %in% methods_to_show]) {
  tables <- get(method)
  cases_linked <- unique(c(tables$all_unsequenced_edges$from, tables$all_unsequenced_edges$to))
  n_cases_linked <- length(cases_linked)
  relevant_cc <- calconnect_dataset %>% 
    filter(CC_RECORD_NUM %in% cases_linked)
  n_sequences_linked <- sum(!is.na(relevant_cc$strain))
  seq_rate <- n_sequences_linked / n_cases_linked
  ci <- binom.test(x = n_sequences_linked, n = n_cases_linked, conf.level = 0.95)$conf.int
  pct_hispanic = mean(relevant_cc$ETHNICITY_SIMPLE == "Hispanic/Latino")
  trow <- tribble(~method, ~num_cases, ~num_sequences, ~seq_rate, ~ci_lower, ~ci_upper, ~pct_hisp,
                  method, n_cases_linked, n_sequences_linked, seq_rate, ci[[1]], ci[[2]], pct_hispanic)
  seq_rate_by_method <- bind_rows(seq_rate_by_method, trow)
}

seq_rate_by_method <- seq_rate_by_method %>% 
  left_join(v4_performance, by=c("method")) %>% 
  arrange(desc(sim)) %>%
  select(
    method,
    num_cases,
    num_sequences,
    seq_rate,
    ci_lower = ci_lower.x,
    ci_upper = ci_upper.x,
    pct_hisp
  ) %>% 
  mutate(method = map_methods(method))
write_csv(seq_rate_by_method, "[REDACTED]")

# Fig S8: Thresholds ----

combined <- bind_rows(
  v4_performance %>% mutate(version="2 Base Pairs"),
  v4_performance_threshold1 %>% mutate(version = "1 Base Pair"),
  v4_performance_threshold3 %>% mutate(version = "3 Base Pairs"),
  v4_performance_threshold4 %>% mutate(version = "4 Base Pairs"),
  v4_performance_threshold5 %>% mutate(version = "5 Base Pairs")
)
version_levels <- rev(c("1 Base Pair","2 Base Pairs","3 Base Pairs", "4 Base Pairs", "5 Base Pairs"))
pdf("[REDACTED]", width=12, height=6)
combined %>% 
  filter(method %in% methods_to_show) %>% 
  ggplot(aes(x = sim, y = factor(map_methods(method), levels=(map_methods(levels))), 
             xmax = ci_upper, xmin=ci_lower, color=factor(version, levels=version_levels))) + 
  geom_point(size=1.5, position=position_dodge2(width=0.9)) +
  geom_errorbar(width=0.6, position=position_dodge(width=0.9)) +
  geom_hline(yintercept = 1:length(methods_to_show) + 0.5, color="gray") +
  scale_color_discrete(guide=guide_legend(reverse=T)) +
  theme_minimal() +
  labs(
    y = "Surveillance Strategy", x = "Informational Value", color="Threshold")
dev.off()

# Fig S9: Alternative Uncertainty Measures ----
combined <- bind_rows(
  v1_performance %>% mutate(uncertainty="Pairwise + Binomial Interval"),
  v2_performance %>% mutate(uncertainty="Pairwise + Case-Bootstrap Interval"),
  v3_performance %>% mutate(uncertainty="Phylo, No Sampling Uncertainty"),
  v4_performance %>% mutate(uncertainty="Phylo + Case-Bootstrap Interval")
)
uncertainty_levels <- rev(c("Phylo, No Sampling Uncertainty", "Pairwise + Binomial Interval",
                            "Pairwise + Case-Bootstrap Interval", "Phylo + Case-Bootstrap Interval"))
pdf("[REDACTED]", width=12, height=6)
combined %>% 
  filter(method %in% methods_to_show) %>% 
  ggplot(aes(x = sim, y = factor(map_methods(method), levels=(map_methods(levels))), 
             xmax = ci_upper, xmin=ci_lower, color=factor(uncertainty, levels=uncertainty_levels))) + 
  geom_point(size=1.5, position=position_dodge2(width=0.9)) +
  geom_errorbar(size=1, width=0.6, position=position_dodge(width=0.9)) +
  geom_hline(yintercept = 1:length(methods_to_show) + 0.5, color="gray") +
  scale_color_discrete(guide=guide_legend(reverse=T)) +
  theme_minimal() +
  labs(
    y = "Surveillance Strategy", x = "Informational Value", color="Methodology")
dev.off()
