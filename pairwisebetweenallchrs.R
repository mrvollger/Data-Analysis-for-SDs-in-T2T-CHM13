


sds = SEDEF[SEDEF$original,]
sds$New = overlap_either(sds, NEW)


p=ggplot(data = sds, aes(color=New, x = aln_len, y=fracMatch)) + 
  geom_point() +
  scale_x_continuous(trans="log10", labels = comma) + 
  #scale_y_continuous(trans="log10", labels = comma) +
  scale_color_manual(values= c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual(values= c(OLDCOLOR, NEWCOLOR))+
  theme_cowplot() +theme(legend.position = "bot")

top = ggplot(data = sds, aes(color=New, fill=New,x = aln_len))+ #y=fracMatch)) + 
  #geom_histogram(aes(x=aln_len), position = "dodge", bins=20) +
  geom_density(alpha = 0.8)+
  scale_x_continuous(trans="log10", labels = comma) + 
  #scale_y_continuous(trans="log10", labels = comma) +
  scale_color_manual(values= c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual(values= c(OLDCOLOR, NEWCOLOR))+
  theme_nothing() + theme(legend.position = "none"); top

right = ggplot(data = sds, aes(color=New, fill=New,y=fracMatch)) + 
  #geom_histogram(aes(x=aln_len), position = "dodge", bins=20) +
  geom_density(alpha = 0.8)+
  #scale_x_continuous(trans="log10", labels = comma) + 
  #scale_y_continuous(trans="log10", labels = comma) +
  scale_color_manual(values= c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual(values= c(OLDCOLOR, NEWCOLOR))+
  theme_nothing() + theme(legend.position = "none"); right

plot_grid(plot_grid(top, NULL, nrow=1, rel_widths = c(4,1)),
          plot_grid(p,  right, nrow=1, rel_widths = c(4,1)),
          rel_heights = c(1,4), ncol=1, align = "v")

  p=ggplot(data=SEDEF)+
  geom_density(aes(x=fracMatch))+
  facet_grid(chr ~ chr2) + theme_cowplot()

ggsave("~/Desktop/tmp.pdf", height = 40, width = 40, plot = p)


p=ggplot(data=SEDEF %>%
           group_by(chr, chr2) %>% summarise(bp_aligned = sum(matchB + mismatchB))
         )+
  geom_bar(aes(y=bp_aligned, group="bp", color = bp_aligned))+
  #scale_x_continuous(trans="log10", labels = comma) + 
  scale_y_continuous(trans="log10", labels = comma) +
  facet_grid(chr ~ chr2) + theme_cowplot(); p

ggsave("~/Desktop/n_mbp.pdf", height = 40, width = 40, plot = p)
