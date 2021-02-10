include("code.jl")
using Measures
pal = (
  C1 = RGB(160/255,47/255,37/255),
  C2 = RGB(230/255,159/255,0/255),
  C3 = RGB(91/255, 135/255, 172/255)
  )


### Fig 01 - PCA with groups of Networks

#### Communities can be clearly grouped according to network betadiversity.
# PCA beta diversity groups

 R"""
   library(ggfortify)
   library(ggplot2)
   library(gridExtra)

   my_pal <- c('black', '#588C7E', '#F2AE72', '#A02F25', '#F2E394', 'ivory2', 'white')
   theme_set(
     theme_minimal() +
       theme(legend.position = "right",
       panel.background = element_blank(),
       panel.grid.minor = element_blank(),
       axis.line = element_line(colour = "ivory2"))
     )

   a <- autoplot(hp_beta_pca, loadings=TRUE, loadings.label = TRUE,
       loadings.label.size=4, loadings.label.colour = 'black',
       loadings.colour = 'black', loadings.label.repel=TRUE,
       main = "(A)") +
       stat_ellipse(geom = "polygon", alpha = 0.2,
       aes(group=kc$cluster, fill = as.factor(kc$cluster)),
         show.legend = FALSE) +
       scale_fill_manual(values=my_pal[2:4]) +
       geom_point(aes(color=as.factor(kc$cluster))) +
       labs(color = "K-means groups") +
       scale_colour_manual(values = my_pal[2:4]) +
       theme(plot.title=element_text(hjust=0))

   # Fleas PCDp
   b <- autoplot(fleas_pcdp_pca, loadings=FALSE,
       main = "(B)") +
       stat_ellipse(geom = "polygon", alpha = 0.2,
       aes(group=kc_fleas_pcdp$cluster, fill = as.factor(kc_fleas_pcdp$cluster)),
         show.legend = FALSE) +
       scale_fill_manual(values=my_pal[2:6]) +
       geom_point(aes(color=as.factor(kc_fleas_pcdp$cluster))) +
       labs(color = "K-means groups") +
       scale_colour_manual(values = my_pal[2:6]) +
       theme(plot.title=element_text(hjust=0))

   # Hosts PCDp
   c <- autoplot(hosts_pcdp_pca, loadings=FALSE,
       main = "(C)") +
       stat_ellipse(geom = "polygon", alpha = 0.2,
       aes(group=kc_hosts_pcdp$cluster, fill = as.factor(kc_hosts_pcdp$cluster)),
         show.legend = FALSE) +
       scale_fill_manual(values=my_pal[2:6]) +
       geom_point(aes(color=as.factor(kc_hosts_pcdp$cluster))) +
       labs(color = "K-means groups") +
       scale_colour_manual(values = my_pal[2:6]) +
       theme(plot.title=element_text(hjust=0))

   # Fleas PCDc
   d <- autoplot(fleas_pcdc_pca, loadings=FALSE,
       main = "(D)") +
       stat_ellipse(geom = "polygon", alpha = 0.2,
       aes(group=kc_fleas_pcdc$cluster, fill = as.factor(kc_fleas_pcdc$cluster)),
         show.legend = FALSE) +
       scale_fill_manual(values=my_pal[2:6]) +
       geom_point(aes(color=as.factor(kc_fleas_pcdc$cluster))) +
       labs(color = "K-means groups") +
       scale_colour_manual(values = my_pal[2:6]) +
       theme(plot.title=element_text(hjust=0))

   # Hosts PCDc
   e <- autoplot(hosts_pcdc_pca, loadings=FALSE,
       main = "(E)") +
       stat_ellipse(geom = "polygon", alpha = 0.2,
       aes(group=kc_hosts_pcdc$cluster, fill = as.factor(kc_hosts_pcdc$cluster)),
         show.legend = FALSE) +
       scale_fill_manual(values=my_pal[2:6]) +
       geom_point(aes(color=as.factor(kc_hosts_pcdc$cluster))) +
       labs(color = "K-means groups") +
       scale_colour_manual(values = my_pal[2:6]) +
       theme(plot.title=element_text(hjust=0))

   png("figs/fig1.png", width = 1024)
   grid.arrange(a, b, c, d, e, layout_matrix = rbind(c(1, 2, 3),
                       c(1, 4, 5)), widths=c(1,1,1))
   dev.off()
 """

### Fig 02 - Panel wih 3 plots beta-diversity vs. PCD
#### Each beta-diversity index relates in a particular way with phylogenetically community dissimilarity (PCD).

# Phylogenetic community dissimilarity vs. first components of each network dissimilarity components analysis

## Species composition dissimilarity vs. first components of each network dissimilarity components analysis, colored by Betadiversity groups
# Fleas PCDc and βs - kmeans of fleas pcdc
a = Plots.scatter(fleas_pcdc_pca.vals[5][:,1], hp_beta_pca_mat[:,1], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, ylabel="\\betas Dissimilarity", guidefontsize=25, size=(790, 526), dpi=300)
# Hosts PCDc and βs - kmeans of hosts pcdc
b = Plots.scatter(hosts_pcdc_pca.vals[5][:,1], hp_beta_pca_mat[:,2], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:bottomright, label = ["\\betawn" "\\betaos" "\\betas"], legendfont=18, size=(790, 556),dpi=300, foreground_color_legend = nothing)

# Fleas PCDc vs. βos
c = Plots.scatter(fleas_pcdc_pca.vals[5][:,1], hp_beta_pca_mat[:,3], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, ylabel="\\betaos Dissimilarity", guidefontsize=25, size=(780, 526))
# Hosts PCDc vs. βos
d = Plots.scatter(hosts_pcdc_pca.vals[5][:,1], hp_beta_pca_mat[:,3], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, size=(780, 526))

# Fleas PCDc vs. βwn
e = Plots.scatter(fleas_pcdc_pca.vals[5][:,1], hp_beta_pca_mat[:,4], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, xlabel="PCDc fleas", ylabel="\\betawn Dissimilarity", guidefontsize=25, size=(780, 526))
# Hosts PCDc vs. βwn
g = Plots.scatter(hosts_pcdc_pca.vals[5][:,1], hp_beta_pca_mat[:,4], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, xlabel="PCDc rodents", guidefontsize=25, size=(780, 526))

#legend = Plots.plot([0 0 0], title = "K-means \n groups", titlefontsize=20, legendfont=18, color = [pal.C1 pal.C2 pal.C3], showaxis = false, grid = false, label = ["\\betawn" "\\betaos" "\\betas"])
title = Plots.plot(title = "Effect of compositional community dissimilarity (PCDc) of parasites (left) and hosts (right) \n on network betadiversity, grouped by kmeans of betadiversity variables", grid = false, showaxis = false, bottom_margin = -150Plots.px, titlefontsize=30)
l = @layout [A{.1h}; Plots.grid(3,2)]

p = Plots.plot(title, a, b, c, d, e, g, size = (2400,2400), layout = l, dpi=300)
display(p)


## Phylogenetic community dissimilarity vs. first components of each network dissimilarity components analysis, colored by Betadiversity groups
# Fleas PCDp and Bs - kmeans of fleas pcdc
a = Plots.scatter(fleas_pcdp_pca.vals[5][:,1], hp_beta_pca_mat[:,1], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, ylabel="\\betas Dissimilarity", guidefontsize=25, dpi=300)
# Hosts PCDp and βs - kmeans of hosts pcdc
b = Plots.scatter(hosts_pcdp_pca.vals[5][:,1], hp_beta_pca_mat[:,2], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:bottomright, label = ["\\betawn" "\\betaos" "\\betas"], legendfont=18, size=(790, 556),dpi=300, foreground_color_legend = nothing)

# Fleas PCDp vs. βos
c = Plots.scatter(fleas_pcdp_pca.vals[5][:,1], hp_beta_pca_mat[:,3], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, ylabel="\\betaos Dissimilarity", guidefontsize=25)
# Hosts PCDp vs. βos
d = Plots.scatter(hosts_pcdp_pca.vals[5][:,1], hp_beta_pca_mat[:,3], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none)

# Fleas PCDp vs. βwn
e = Plots.scatter(fleas_pcdp_pca.vals[5][:,1], hp_beta_pca_mat[:,4], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, xlabel="PCDp fleas", ylabel="\\betawn Dissimilarity", guidefontsize=25)
# Hosts PCDp vs. βwn
g = Plots.scatter(hosts_pcdp_pca.vals[5][:,1], hp_beta_pca_mat[:,4], groups=kc.vals[1], color = [pal.C1 pal.C2 pal.C3], ms=5, markerstrokewidth=0, xlim = [-8,9], ylim = [-2.5, 2.5], margin=20mm, legend=:none, xlabel="PCDp rodents", guidefontsize=25)

#legend = Plots.plot([0 0 0], title = "K-means \n groups", titlefontsize=20, legendfont=18, color = [pal.C1 pal.C2 pal.C3], showaxis = false, grid = false, label = ["\\betawn" "\\betaos" "\\betas"])
title = Plots.plot(title = "Effect of phylogenetic community dissimilarity (PCDp) of parasites (left) and hosts (right) \n on network betadiversity, grouped by kmeans of betadiversity variables", grid = false, showaxis = false, bottom_margin = -150Plots.px, titlefontsize=30)
l = @layout [A{.1h}; Plots.grid(3,2)]

p = Plots.plot(title, a, b, c, d, e, g, size = (2400,2400), layout = l)
display(p)

### Fig 03 - Map
#### The separation of communities by components of beta-diversity was also observed geographically
beta_map = [coordinates kc.vals[1] kc_fleas_pcdp.vals[1] kc_hosts_pcdp.vals[1] kc_fleas_pcdc.vals[1] kc_hosts_pcdc.vals[1]];
range_map = [(minimum(coordinates[:,1])-15, maximum(coordinates[:,1])+15), (minimum(coordinates[:,2])-5, maximum(coordinates[:,2])+15)];
using SimpleSDMLayers
bioclim = SimpleSDMLayers.worldclim(1, resolution = 5.0);
bioclim_copy= Base.similar(bioclim)
filter(!isnothing, bioclim_copy.grid) .= 0

# Betadiversity - Fig 3A
function f(x,n)
    y = x[x[:,3].==n,:]
    return y[:,1], y[:,2]
end

p = heatmap(bioclim_copy[range_map[1], range_map[2]], c=:amp, aspectratio=92.60/60.75)
scatter!(p, f(beta_map, 1), color = [pal.C1], ms=5, label = "", title = "\\betawn", markerstrokewidth=0, dpi=300)
q = heatmap(bioclim_copy[range_map[1], range_map[2]], c=:amp, aspectratio=92.60/60.75)
scatter!(q, f(beta_map, 2), color = [pal.C2], ms=5, label = "", title = "\\betaos", markerstrokewidth=0, dpi=300)
r = heatmap(bioclim_copy[range_map[1], range_map[2]], c=:amp, aspectratio=92.60/60.75)
scatter!(r, f(beta_map, 3), color = [pal.C3], ms=5, label = "", title = "\\betas", markerstrokewidth=0, dpi=300)

title = plot(title = "(A) Spatial clusters of networks by beta-diversity metrics", grid = false, showaxis = false, bottom_margin = -50Plots.px)
l = @layout [A{.1h}; grid(3,1)]

plot(title, p, q, r, size = (1024,824),layout = l)


# PCD maps - Fig 3B
#---PCDp MAPS ---
function f(x, g, n)
    y = x[x[:,g].==n,:]
    return y[:,1], y[:,2]
end
using StatsPlots

# Group 1
h_plot_series = [f(beta_map, 4, 1)[1] f(beta_map, 4, 1)[2]; f(beta_map, 5, 3)[1] f(beta_map, 5, 3)[2]]

h = marginalkde(h_plot_series[:,1], h_plot_series[:,2], size = (980,520), colour = "black")
h.subplots[2].series_list[1][:linecolor] = RGBA{Float64}(1.0,1.0,1.0)
xaxis!(h, range_map[1], subplot=1)
yaxis!(h, range_map[2], subplot=3)
plot!(h, bioclim_copy[left=range_map[1][1],right=range_map[1][2], bottom=range_map[2][1], top=range_map[2][2]], c=:amp, subplot=2, aspectratio=:auto)
plot!(h, showaxis=false, subplot=1)
plot!(h, showaxis=false, subplot=3)

scatter!(h, f(beta_map, 4, 1), color = [pal.C3], ms=5, label = "", markerstrokewidth=0, subplot=2) #fleas
scatter!(h, f(beta_map, 5, 3), color = [pal.C2], marker=:xcross, ms=5, label = "", markerstrokewidth=0, dpi=300, subplot=2) #mammals

# Group 2
t_plot_series = [f(beta_map, 4, 2)[1] f(beta_map, 4, 2)[2]; f(beta_map, 5, 1)[1] f(beta_map, 5, 1)[2]]

t = marginalkde(t_plot_series[:,1], t_plot_series[:,2], size = (980,520), colour = "black")
t.subplots[2].series_list[1][:linecolor] = RGBA{Float64}(1.0,1.0,1.0)
xaxis!(t, range_map[1], subplot=1)
yaxis!(t, range_map[2], subplot=3)
plot!(t, bioclim_copy[left=range_map[1][1],right=range_map[1][2], bottom=range_map[2][1], top=range_map[2][2]], c=:amp, subplot=2, aspectratio=:auto)
plot!(t, showaxis=false, subplot=1, title = "PCDp")
plot!(t, showaxis=false, subplot=3)

scatter!(t, f(beta_map, 4, 2), color = [pal.C3], ms=5, label = "", markerstrokewidth=0, subplot=2) #fleas
scatter!(t, f(beta_map, 5, 1), color = [pal.C2], marker=:xcross, ms=5, label = "", markerstrokewidth=0, dpi=300, subplot=2) #mammals

# Group 3
u_plot_series = [f(beta_map, 4, 3)[1] f(beta_map, 4, 3)[2]; f(beta_map, 5, 2)[1] f(beta_map, 5, 2)[2]]

u = marginalkde(u_plot_series[:,1], u_plot_series[:,2], size = (980,520), colour = "black")
u.subplots[2].series_list[1][:linecolor] = RGBA{Float64}(1.0,1.0,1.0)
xaxis!(u, range_map[1], subplot=1)
yaxis!(u, range_map[2], subplot=3)
plot!(u, bioclim_copy[left=range_map[1][1],right=range_map[1][2], bottom=range_map[2][1], top=range_map[2][2]], c=:amp, subplot=2, aspectratio=:auto)
plot!(u, showaxis=false, subplot=1)
plot!(u, showaxis=false, subplot=3)

scatter!(u, f(beta_map, 4, 3), color = [pal.C3], ms=5, label = "", markerstrokewidth=0, subplot=2) #fleas
scatter!(u, f(beta_map, 5, 2), color = [pal.C2], marker=:xcross, ms=5, label = "", markerstrokewidth=0, dpi=300, subplot=2) #mammals


#---PCDc MAPS ---
# Group 1
y_plot_series = [f(beta_map, 6, 1)[1] f(beta_map, 6, 1)[2]; f(beta_map, 7, 3)[1] f(beta_map, 7, 3)[2]]

y = marginalkde(y_plot_series[:,1], y_plot_series[:,2], size = (980,520), colour = "black")
y.subplots[2].series_list[1][:linecolor] = RGBA{Float64}(1.0,1.0,1.0)
xaxis!(y, range_map[1], subplot=1)
yaxis!(y, range_map[2], subplot=3)
plot!(y, bioclim_copy[left=range_map[1][1],right=range_map[1][2], bottom=range_map[2][1], top=range_map[2][2]], c=:amp, subplot=2, aspectratio=:auto)
plot!(y, showaxis=false, subplot=1)
plot!(y, showaxis=false, subplot=3)

scatter!(y, f(beta_map, 6, 1), color = [pal.C3], ms=5, label = "", markerstrokewidth=0, subplot=2) #fleas
scatter!(y, f(beta_map, 7, 3), color = [pal.C2], marker=:xcross, ms=5, label = "", markerstrokewidth=0, dpi=300, subplot=2) #mammals

#Group 2
z_plot_series = [f(beta_map, 6, 2)[1] f(beta_map, 6, 2)[2]; f(beta_map, 7, 1)[1] f(beta_map, 7, 1)[2]]

z = marginalkde(z_plot_series[:,1], z_plot_series[:,2], size = (980,520), colour = "black")
z.subplots[2].series_list[1][:linecolor] = RGBA{Float64}(1.0,1.0,1.0)
xaxis!(z, range_map[1], subplot=1)
yaxis!(z, range_map[2], subplot=3)
plot!(z, bioclim_copy[left=range_map[1][1],right=range_map[1][2], bottom=range_map[2][1], top=range_map[2][2]], c=:amp, subplot=2, aspectratio=:auto)
plot!(z, showaxis=false, subplot=1, title = "PCDc")
plot!(z, showaxis=false, subplot=3)

scatter!(z, f(beta_map, 6, 2), color = [pal.C3], ms=5, label = "", markerstrokewidth=0, subplot=2) #fleas
scatter!(z, f(beta_map, 7, 1), color = [pal.C2], marker=:xcross, ms=5, label = "", markerstrokewidth=0, dpi=300, subplot=2) #mammals


#Group 3
a_plot_series = [f(beta_map, 6, 3)[1] f(beta_map, 6, 3)[2]; f(beta_map, 7, 2)[1] f(beta_map, 7, 2)[2]]

a = marginalkde(a_plot_series[:,1], a_plot_series[:,2], size = (980,520), colour = "black")
a.subplots[2].series_list[1][:linecolor] = RGBA{Float64}(1.0,1.0,1.0)
xaxis!(a, range_map[1], subplot=1)
yaxis!(a, range_map[2], subplot=3)
plot!(a, bioclim_copy[left=range_map[1][1],right=range_map[1][2], bottom=range_map[2][1], top=range_map[2][2]], c=:amp, subplot=2, aspectratio=:auto)
plot!(a, showaxis=false, subplot=1)
plot!(a, showaxis=false, subplot=3)

scatter!(a, f(beta_map, 6, 3), color = [pal.C3], ms=5, label = "", markerstrokewidth=0, subplot=2) #fleas
scatter!(a, f(beta_map, 7, 2), color = [pal.C2], marker=:xcross, ms=5, label = "", dpi=300, subplot=2) #mammals

# final plot
legend = Plots.scatter([0 0], marker=[:circle :xcross], markerstrokewidth=0, legendfont=18, color = [pal.C3 pal.C2], foreground_color_legend = nothing, showaxis = false, grid = false, label = ["Parasites" "Hosts"])
title = plot(title = "Spatial distribution of PCD components clusters for hosts and parasites", grid = false, showaxis = false, bottom_margin = -50Plots.px, titlefontsize=30)
l = @layout [A{.1h}; grid(3,2) E{.15w}]

plot(title, z, t, a, u, y, h, legend, size = (2400,2000),layout = l)
