
# write table for PerBase methylation in SDs
write.table(chm13.ovl, file = paste0(dat, "/sd/t2t_chm13v1.0_SD_clustered_methylationPerBase.bed"),quote = F, sep = "\t", row.names = F, col.names = T)
# make distance labels
chm13.ovl.labs <- chm13.ovl %>%
  group_by(dist) %>%
  summarise(med_meth = median(methylated_frequency)) %>%
  distinct()
# cluster and make distance matrix
clust_meth <- chm13.ovl %>% 
  select(c(ID,methylated_frequency,dist)) %>%
  group_by(ID, dist) %>%
  summarise(meth = mean(methylated_frequency)) %>%
  spread(dist,meth) %>%
  column_to_rownames("ID") %>%
  as.matrix() 
# annotation of flank vs SD for plot
annot <- chm13.ovl %>% 
  select(c(ID,methylated_frequency,dist)) %>%
  group_by(ID, dist) %>%
  summarise(meth = mean(methylated_frequency))%>%
  mutate(reg = case_when( dist < 0 ~ "flank", 
                          dist > bodylen ~ "flank", 
                          TRUE ~ "SD")) %>%
  ungroup() %>%
  select(c(dist, reg)) %>%
  distinct() %>%
  column_to_rownames("dist")
```

```{r, echo=F, fig.height=5, fig.width=8}
# make heatmap and clustered bed
library(viridis)
res <- pheatmap(clust_meth, cluster_cols = F, annotation_col = annot, show_rownames = F, show_colnames = F, main = "SD Methylation")
res
SDmeth.clust <- as.data.frame(clust_meth) %>%
  mutate(avgmeth = rowMeans(., na.rm = T)) %>%
  mutate(cluster_num = cutree(res$tree_row, k = 2)) %>%
  rownames_to_column("ID") %>%
  select(c(ID, cluster_num, avgmeth)) %>%
  mutate(clust_meth = case_when(cluster_num == 1 ~ "unmeth", 
                                cluster_num == 2 ~ "meth"))

final_clust <- merge(flank,SDmeth.clust, by = "ID") %>%
  select(-c(ID)) %>%
  rename(seqnames = "chr")
write.table(final_clust, file = paste0(dat, "/sd/t2t_chm13v1.0_SD_clustered_methylation.bed"),quote = F, sep = "\t", row.names = F, col.names = T)
ggsave(
  paste0(figs, "/","SD_methylation_heatmap.pdf"),
  plot = res,
  scale = 1,
  width = 10,
  height = 10,
)



