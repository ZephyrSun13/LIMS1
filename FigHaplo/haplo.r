
library(survival)
library(survminer)
library(RColorBrewer)
library(ggpubr)


Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

Clinical <- read.delim(file = "Clinical_Haplo.xls", row.names = 1, as.is = T)

Clinical$GL_time <- Clinical$GL_time/30

Clinical$GenoGroupNEJHaplo <- Clinical$GenoGroup

Clinical$GenoGroupNEJHaplo[which(Clinical$NEJGenoGroup == "Risk")] <- "Risk"

Clinical$GenoGroupNEJHaplo[which(is.na(Clinical$NEJGenoGroup))] <- NA

Clinical$GenoGroupNEJHaplo <- factor(Clinical$GenoGroupNEJHaplo, levels = c("Maj", "Min", "Risk", "Other"))


Clinical.CC <- Clinical[Clinical$R2R == "Caucasian > Caucasian", ]

Clinical.noCC <- Clinical[Clinical$R2R != "Caucasian > Caucasian", ]


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroup, data = Clinical)

pp <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
	  pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Red", "Blue"), ],
          legend.title = "Allele type introduced",
          legend.labs = c("Major", "Minor", "Other"),
	  xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
	  font.title = c(10, "bold", "darkblue"), 
          risk.table.fontsize = 4, 
	  risk.table.y.text = F
)

pp$plot <- pp$plot + theme(plot.margin = margin(0, 5, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp$table <- pp$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroup, data = Clinical.CC)

pp2 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Red", "Blue"), ],
          legend.title = "Allele type introduced\n(E -> E)",
          legend.labs = c("Major", "Minor", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp2$plot <- pp2$plot + theme(plot.margin = margin(0, 25, 0, 18, "pt")) + guides(colour = guide_legend(nrow = 2))

pp2$table <- pp2$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ DR_Haplo1_MM_Cat, data = Clinical)

pp3 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("HeatRd3", "HeatRd5", "HeatRd7"), ],
          legend.title = "Number of mismatches",
          legend.labs = c("0", "1", "2"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp3$plot <- pp3$plot + theme(plot.margin = margin(0, 25, 0, 18, "pt")) + guides(colour = guide_legend(nrow = 2))

pp3$table <- pp3$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ DR_Haplo1_MM_Cat, data = Clinical.CC)

pp4 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("HeatRd3", "HeatRd5", "HeatRd7"), ],
          legend.title = "Number of mismatches\n(E -> E)",
          legend.labs = c("0", "1", "2"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp4$plot <- pp4$plot + theme(plot.margin = margin(0, 25, 0, 18, "pt")) + guides(colour = guide_legend(nrow = 2))

pp4$table <- pp4$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical)

pp5 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Red", "SteelDark", "Blue"), ],
          legend.title = "Allele type introduced\n(All)",
          legend.labs = c("Major", "Minor", "rs893403", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp5$plot <- pp5$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp5$table <- pp5$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical[Clinical$GenoGroupNEJHaplo %in% c("Other", "Maj"), ])

pp5.1 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Blue"), ],
          legend.title = "Allele type introduced\n(All)",
          legend.labs = c("Major", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp5.1$plot <- pp5.1$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp5.1$table <- pp5.1$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical[Clinical$GenoGroupNEJHaplo %in% c("Other", "Min"), ])

pp5.2 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Red", "Blue"), ],
          legend.title = "Allele type introduced\n(All)",
          legend.labs = c("Minor", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp5.2$plot <- pp5.2$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp5.2$table <- pp5.2$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical[Clinical$GenoGroupNEJHaplo %in% c("Other", "Risk"), ])

pp5.3 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("SteelDark", "Blue"), ],
          legend.title = "Allele type introduced\n(All)",
          legend.labs = c("rs893403", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp5.3$plot <- pp5.3$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp5.3$table <- pp5.3$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.CC)

pp6 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Red", "SteelDark", "Blue"), ],
          legend.title = "Allele type introduced\n(E -> E)",
          legend.labs = c("Major", "Minor", "rs893403", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp6$plot <- pp6$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp6$table <- pp6$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.CC[Clinical.CC$GenoGroupNEJHaplo %in% c("Maj", "Other"), ])

pp6.1 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Blue"), ],
          legend.title = "Allele type introduced\n(E -> E)",
          legend.labs = c("Major", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp6.1$plot <- pp6.1$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp6.1$table <- pp6.1$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.CC[Clinical.CC$GenoGroupNEJHaplo %in% c("Min", "Other"), ])

pp6.2 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Red", "Blue"), ],
          legend.title = "Allele type introduced\n(E -> E)",
          legend.labs = c("Minor", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp6.2$plot <- pp6.2$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp6.2$table <- pp6.2$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.CC[Clinical.CC$GenoGroupNEJHaplo %in% c("Risk", "Other"), ])

pp6.3 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("SteelDark", "Blue"), ],
          legend.title = "Allele type introduced\n(E -> E)",
          legend.labs = c("rs893403", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp6.3$plot <- pp6.3$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp6.3$table <- pp6.3$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")

#ggsave("FigHaplo.pdf", ggarrange(ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), ggarrange(pp4$plot, pp4$table, nrow = 2, heights = c(2, 1)), ggarrange(pp$plot, pp$table, nrow = 2, heights = c(2, 1)), ggarrange(pp2$plot, pp2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp5$plot, pp5$table, nrow = 2, heights = c(2, 1)), ggarrange(pp6$plot, pp6$table, nrow = 2, heights = c(2, 1)), ncol = 2, nrow = 3), height = 15, width = 8)

#ggsave("FigHaplo.pdf", ggarrange(ggarrange(ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), ggarrange(pp4$plot, pp4$table, nrow = 2, heights = c(2, 1)), nrow = 2), ggarrange(ggarrange(pp$plot, pp$table, nrow = 2, heights = c(2, 1)), ggarrange(pp2$plot, pp2$table, nrow = 2, heights = c(2, 1)), nrow = 2), ggarrange(ggarrange(pp5$plot, pp5$table, nrow = 2, heights = c(2, 1)), ggarrange(pp6$plot, pp6$table, nrow = 2, heights = c(2, 1)), nrow = 2), labels = c("A", "B", "C"), ncol = 3, nrow = 1), height = 10, width = 12)

#ggsave("FigHaplo.pdf", ggarrange(ggarrange(ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), ggarrange(pp4$plot, pp4$table, nrow = 2, heights = c(2, 1)), nrow = 2, labels = c("A", "B")), ggarrange(ggarrange(pp5$plot, pp5$table, nrow = 2, heights = c(2, 1)), ggarrange(pp6$plot, pp6$table, nrow = 2, heights = c(2, 1)), nrow = 2, labels = c("C", "D")), ncol = 2, nrow = 1), height = 10, width = 8)

ggsave("FigHaplo.pdf", ggarrange(ggarrange(ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), ggarrange(pp4$plot, pp4$table, nrow = 2, heights = c(2, 1)), nrow = 2), ggarrange(ggarrange(pp5$plot, pp5$table, nrow = 2, heights = c(2, 1)), ggarrange(pp6$plot, pp6$table, nrow = 2, heights = c(2, 1)), nrow = 2), ncol = 2, nrow = 1), height = 10, width = 8)


ggsave("FigHaplo.panel1.pdf", ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), height = 5, width = 4)

ggsave("FigHaplo.panel4.pdf", ggarrange(pp6$plot, pp6$table, nrow = 2, heights = c(2, 1)), height = 5, width = 4)


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.noCC)

pp7 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Red", "Blue", "SteelDark"), ],
          legend.title = "Allele type introduced\n(non E -> E)",
          legend.labs = c("Major", "Minor", "rs893403", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp7$plot <- pp7$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp7$table <- pp7$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")

ggsave("FigHaploNoCC.pdf", ggarrange(pp7$plot, pp7$table, nrow = 2, heights = c(2, 1)), height = 5, width = 4)


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.noCC[Clinical.noCC$GenoGroupNEJHaplo %in% c("Other", "Maj"), ])

pp7.1 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Green", "Blue"), ],
          legend.title = "Allele type introduced\n(non E -> E)",
          legend.labs = c("Major", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp7.1$plot <- pp7.1$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp7.1$table <- pp7.1$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.noCC[Clinical.noCC$GenoGroupNEJHaplo %in% c("Other", "Min"), ])

pp7.2 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("Red", "Blue"), ],
          legend.title = "Allele type introduced\n(non E -> E)",
          legend.labs = c("Minor", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp7.2$plot <- pp7.2$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp7.2$table <- pp7.2$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


fit <- survfit(Surv(GL_time, GL_event) ~ GenoGroupNEJHaplo, data = Clinical.noCC[Clinical.noCC$GenoGroupNEJHaplo %in% c("Other", "Risk"), ])

pp7.3 <- ggsurvplot(fit,
          pval = TRUE, conf.int = F,
          pval.size = 4,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_classic(), # Change ggplot2 theme
          palette = Colors[c("SteelDark", "Blue"), ],
          legend.title = "Allele type introduced\n(non E -> E)",
          legend.labs = c("rs893403", "Other"),
          xlab = "Follow-up time (months)",
          fontsize = 6,
          font.x = 10,
          font.y = 10,
          font.legend = 10,
          font.tickslab = 10,
          font.title = c(10, "bold", "darkblue"),
          risk.table.fontsize = 4,
          risk.table.y.text = F
)

pp7.3$plot <- pp7.3$plot + theme(plot.margin = margin(0, 25, 0, 18)) + guides(colour = guide_legend(nrow = 2))

pp7.3$table <- pp7.3$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")



ggsave("FigHaplo.sep.pdf", ggarrange(ggarrange(ggarrange(pp5.1$plot, pp5.1$table, nrow = 2, heights = c(2, 1)), ggarrange(pp5.2$plot, pp5.2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp5.3$plot, pp5.3$table, nrow = 2, heights = c(2, 1)), ncol = 3, labels = c("A", "B", "C")), ggarrange(ggarrange(pp6.1$plot, pp6.1$table, nrow = 2, heights = c(2, 1)), ggarrange(pp6.2$plot, pp6.2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp6.3$plot, pp6.3$table, nrow = 2, heights = c(2, 1)), ncol = 3, labels = c("D", "E", "F")), ggarrange(ggarrange(pp7.1$plot, pp7.1$table, nrow = 2, heights = c(2, 1)), ggarrange(pp7.2$plot, pp7.2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp7.3$plot, pp7.3$table, nrow = 2, heights = c(2, 1)), ncol = 3, labels = c("G", "H", "I")), ncol = 1, nrow = 3), height = 15, width = 12)


ggsave("FigHaplo.nonEE.sep.pdf", ggarrange(ggarrange(pp7.1$plot, pp7.1$table, nrow = 2, heights = c(2, 1)), ggarrange(pp7.2$plot, pp7.2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp7.3$plot, pp7.3$table, nrow = 2, heights = c(2, 1)), ncol = 3, labels = c("A", "B", "C")), height = 5, width = 12)


	Clinical$GenoGroup <- factor(Clinical$GenoGroup, levels = c("Other", "Min", "Maj"))

        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "GenoGroup", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum + NEJRiskHeta")), data = Clinical)

	Result <- cbind(paste0("(", round(summary(fit)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit)$coefficients[, c(2, 5)], 2), summary(fit)$n) 

        fit.cc <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "GenoGroup", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum + NEJRiskHeta")), data = Clinical[Clinical$R2R == "Caucasian > Caucasian", ])

        Result.cc <- cbind(paste0("(", round(summary(fit.cc)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit.cc)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit.cc)$coefficients[, c(2, 5)], 2), summary(fit.cc)$n)

	write.table(rbind(Result, Result.cc), file = "GenoGroup.xls", sep = "\t", quote = F, col.names = NA)


        Clinical$GenoGroupNEJHaplo <- factor(Clinical$GenoGroupNEJHaplo, levels = c("Other", "Risk", "Min", "Maj"))

        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "GenoGroupNEJHaplo", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum")), data = Clinical)

        Result <- cbind(paste0("(", round(summary(fit)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit)$coefficients[, c(2, 5)], 2), summary(fit)$n)

        fit.cc <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "GenoGroupNEJHaplo", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum")), data = Clinical[Clinical$R2R == "Caucasian > Caucasian", ])

        Result.cc <- cbind(paste0("(", round(summary(fit.cc)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit.cc)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit.cc)$coefficients[, c(2, 5)], 2), summary(fit.cc)$n)

        fit.nocc <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "GenoGroupNEJHaplo", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum")), data = Clinical[Clinical$R2R != "Caucasian > Caucasian", ])

        Result.nocc <- cbind(paste0("(", round(summary(fit.nocc)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit.nocc)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit.nocc)$coefficients[, c(2, 5)], 2), summary(fit.nocc)$n)

        write.table(rbind(Result, Result.cc, Result.nocc), file = "GenoGroupNEJHaplo.xls", sep = "\t", quote = F, col.names = NA)


        Clinical$DR_Haplo1_MM_Cat <- factor(Clinical$DR_Haplo1_MM_Cat, levels = c("0", "1", "2"))

        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "DR_Haplo1_MM_Cat", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum + NEJRiskHeta")), data = Clinical)

        Result <- cbind(paste0("(", round(summary(fit)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit)$coefficients[, c(2, 5)], 2), summary(fit)$n)

        fit.cc <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "DR_Haplo1_MM_Cat", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum + NEJRiskHeta")), data = Clinical[Clinical$R2R == "Caucasian > Caucasian", ])

        Result.cc <- cbind(paste0("(", round(summary(fit.cc)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit.cc)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit.cc)$coefficients[, c(2, 5)], 2), summary(fit.cc)$n)

        write.table(rbind(Result, Result.cc), file = "DR_Haplo1_MM_Cat.xls", sep = "\t", quote = F, col.names = NA)

