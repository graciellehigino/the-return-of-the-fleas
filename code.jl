using Pkg
Pkg.activate(".")
Pkg.instantiate

using EcologicalNetworks
using JLD
using LinearAlgebra
using Mangal
using Plots
using Statistics


hp_data = dataset("hadfield_2014")
hp_networks = networks(hp_data)

hp_taxa_array = vcat(nodes.(hp_networks)...);
N = [convert(Mangal.UnipartiteNetwork, n) for n in hp_networks];
Nt = taxonize.(N);
hp_taxa_array = Any[]
for i in 1:length(Nt)
    hp_taxa_subarray = String[]
    for j in Nt[i]
      push!(hp_taxa_subarray, "$(j.from.name)")
      push!(hp_taxa_subarray, "$(j.to.name)")
    end
    push!(hp_taxa_array, hp_taxa_subarray)
end

## Phylogenies
# How to get the phylogenies from DataDryad?
# using ZipFile
# using HTTP
# using Phylo
using RCall

@rput hp_taxa_array

begin
    R"""
    library(ape)
    tree_fleas <- read.tree("data/FleaTree.tre")
    tree_hosts <- read.tree("data/HostTree.tre")

    # Communities matrices
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(reshape2)
    library(vegan)

    hp_df <- map(map(hp_taxa_array, unique), strsplit, "\t")
    hp_df <- melt(hp_df)
    hp_df <- pivot_wider(hp_df, names_from=value, values_from=L2)
    hp_sites <- hp_df[,1]
    hp_df <- select(hp_df, -1)
    hp_df[!is.na(hp_df)]<-1
    hp_df[is.na(hp_df)]<-0
    #hp_df[,1] <- hp_sites
    hp_mat<-as.matrix(hp_df)
    colnames(hp_mat) <- gsub(" ", "_", colnames(hp_mat))

    #Inspect sites with no representatives in phylogeny
    # hosts
    hp_mat2 = hp_mat[, colnames(hp_mat) %in% tree_hosts$tip.label]
    drop_hosts = which(rowSums(hp_mat2) == 0)

    # parasites
    hp_mat2 = hp_mat[, colnames(hp_mat) %in% tree_fleas$tip.label]
    drop_parasites = which(rowSums(hp_mat2) == 0)

    hp_mat2 = hp_mat[-c(drop_hosts, drop_parasites),]

    # PD of hosts and parasites
    library(phyr)

    # PD of communities
    pcd_hosts <- pcd(hp_mat2, tree_hosts)
    pcd_fleas <- pcd(hp_mat2, tree_fleas)
    hosts_pcdp <- pcd_hosts$PCDp
    hosts_pcdc <- pcd_hosts$PCDc
    fleas_pcdp <- pcd_fleas$PCDp
    fleas_pcdc <- pcd_fleas$PCDc

    # Summarizing phylogenetic community diversity
    hosts_pcdp <- replace(hosts_pcdp, is.infinite(hosts_pcdp),0) # replace Inf by 0; probably the site where there's no correspondence in tree hosts
    hosts_pcdp_pca <- prcomp(hosts_pcdp, scale = TRUE)
    hosts_pcdc_pca <- prcomp(hosts_pcdc, scale = TRUE)
    fleas_pcdp_pca <- prcomp(fleas_pcdp, scale = TRUE)
    fleas_pcdc_pca <- prcomp(fleas_pcdc, scale = TRUE)
    """
end
@rget drop_hosts
@rget drop_parasites
@rget fleas_pcdp
@rget fleas_pcdc
@rget hosts_pcdp
@rget hosts_pcdc
@rget hosts_pcdp_pca
@rget hosts_pcdc_pca
@rget fleas_pcdp_pca
@rget fleas_pcdc_pca

# PD betadiversity
using PhyloNetworks, Diversity

deleteat!(Nt, drop_hosts)
deleteat!(Nt, drop_parasites)
# Inspect the metanetwork
M = reduce(∪, Nt)

# connectance
connectance(M)

# modularity


# nestedness
B = [convert(BipartiteNetwork, n) for n in Nt];
η(reduce(∪, B))

# Inspect realisations
# For the following metrics, a = species or links in common between networks,
# b = a - species or links from j, and c = a - species or links from i.

Bs = KGL02.([EcologicalNetworks.βs(i, j) for i in Nt, j in Nt]);      # dissimilarity in the species composition of communities

# Bs for hosts or parasites
function split_βs(X::T, Y::T; dims::Int64) where {T<:BinaryNetwork}
    a = richness(intersect(X,Y); dims = dims)
    b = richness(Y; dims = dims)-a
    c = richness(X; dims = dims)-a
    return (a=a, b=b, c=c)
end

# Hosts
Bs_hosts = KGL11.([split_βs(i, j; dims = 1) for i in B, j in B]);

# Parasites
Bs_fleas = KGL11.([split_βs(i, j; dims = 2) for i in B, j in B]);

Bos = replace!(KGL11.([EcologicalNetworks.βos(i, j) for i in Nt, j in Nt]), NaN => 0);    # dissimilarity of interactions in co-occuring species
Bwn = KGL11.([EcologicalNetworks.βwn(i, j) for i in Nt, j in Nt]);    # differences between interactions
Bst = Bwn.-Bos;
BosT = KGL02.([βos(M, i) for i in Nt]);
# C_local = connectance.(Nt);
# Nest_local = η.(B);

## Summarizing B-diversity variables

# PCoA
hp_beta_list = [Bs_fleas, Bs_hosts, Bos, Bwn]

@rput hp_beta_list
begin
    R"""
    hp_beta_pca_mat <- matrix(0, nrow=nrow(hp_mat2),length(hp_beta_list))
    for (i in 1:length(hp_beta_list)){
        hp_beta_pca_mat[,i] <- prcomp(hp_beta_list[[i]])$x[,1]
        colnames(hp_beta_pca_mat) <- c("βs_parasites", "βs_hosts", "βos", "βwn")
    }
    hp_beta_pca <- prcomp(hp_beta_pca_mat, scale = TRUE)

    set.seed(42)
    kc <- kmeans(hp_beta_pca$x, 3, iter=50)
    kc_hosts_pcdp <- kmeans(hosts_pcdp_pca$x, 3, iter=50)
    kc_hosts_pcdc <- kmeans(hosts_pcdc_pca$x, 3, iter=50)
    kc_fleas_pcdp <- kmeans(fleas_pcdp_pca$x, 3, iter=50)
    kc_fleas_pcdc <- kmeans(fleas_pcdc_pca$x, 3, iter=50)
    """
end
@rget hp_beta_pca
@rget hp_beta_pca_mat
@rget kc
@rget kc_fleas_pcdc
@rget kc_fleas_pcdp
@rget kc_hosts_pcdp
@rget kc_hosts_pcdc


# Geographic distance
using Distances

coordinates = zeros(length(hp_networks),2);
for i in 1:length(hp_networks)
    coordinates[i,1] = hp_networks[i].position.coordinates[1]
    coordinates[i,2] = hp_networks[i].position.coordinates[2]
end
coordinates = coordinates[1:end .!=22, 1:end];
