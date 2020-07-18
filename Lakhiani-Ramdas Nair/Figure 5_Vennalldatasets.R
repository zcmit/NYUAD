source('code/Functions.R')

# Figure 6 - Venns --------------------------------------------------------

fig6_targets_venn = make_3way_eulerr(intersect(targets,natf6_sig$ensembl_id),
                                     intersect(targets,ias_sig$ensembl_id),
                                     intersect(targets,etoh_sig$ensembl_id),
                                     'figures/Fig 6/fig6_targets_venn',
                                     c(natf6_colour,ias_colour,etoh_colour),7)

fig6_upr_venn = make_3way_eulerr(intersect(upr,natf6_sig$ensembl_id),
                                     intersect(upr,ias_sig$ensembl_id),
                                     intersect(upr,etoh_sig$ensembl_id),
                                     'figures/Fig 6/fig6_upr_venn',
                                     c(natf6_colour,ias_colour,etoh_colour),7)
