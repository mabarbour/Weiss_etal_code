
## Setup ----

# load libraries
library(tidyverse)
library(cowplot)

# ggplot theme
theme_tess <- function () { 
  theme_cowplot()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))+
    theme(axis.text.x=element_text(size=20))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.x=element_text(size=20))+
    theme(axis.title.y=element_text(size=20))+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    theme(axis.title.y=element_text(size=20))
}

data <-read.csv("data/body_size_final.csv",stringsAsFactors = FALSE,
                strip.white = TRUE, na.strings = c("NA","") )

summary_data <- data %>%
  group_by(temperature,flour,population_replicate) %>%
  group_by(temperature,flour)%>%
  summarise(
    mean = mean(weightingrams,na.rm = TRUE),
    se = sd(weightingrams, na.rm = TRUE) / sqrt(sum(!is.na(weightingrams)))
  )

B_G_VG_pal <- c("black","orange", "lightgreen", "darkgreen")

p<-ggplot(summary_data, aes(x = temperature, y = mean,color=flour)) +
  geom_point(size = 3,position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.3),
                width = 0) +
  xlab("Temperature (\u00B0C)") +
  ylab("Weight (g)")+
  scale_x_discrete(limits = c("founder", "27.5", "30","32.5","35"),
                   labels=c("Founder", "27.5", "30","32.5","35"))+
  scale_color_manual(values = c("black","orange", "lightgreen", "darkgreen"), 
                     name = "Resource quality", 
                     breaks=c("founder", "B", "G", "VG"),
                     labels = c("Founder","Bad","Medium","Good"),
                     guide = guide_legend(reverse=TRUE))+
  theme_tess()
  #annotate("text",label =expression(bold("Temperature " ~bolditalic("P")~" < 0.001")), 
           #x=4,y=0.00109,size=4.6)+
  #annotate("text",label =expression(bold("Resource NS")), x=4,y=0.00108,size=4.6)+
  #annotate("text",label =expression(bold("Temperature*Resource NS")), x=4,y=0.00107,size=4.6)+

p

ggsave(file="figures/Body size.pdf", p, width = 22, height = 21, units = "cm")

#### ANALYSIS #####

#need to get means of all the individuals measured in each population first
summary_data_2 <- data %>%
  filter(flour!="VB")%>%
  group_by(temperature,flour,population_replicate) %>%
  summarise(
    mean = mean(weightingrams,na.rm = TRUE),
    se = sd(weightingrams, na.rm = TRUE) / sqrt(sum(!is.na(weightingrams)))
  )

summary_data_2<-subset (summary_data_2,temperature!="founder")
summary_data_2$temperature<-as.numeric(summary_data_2$temperature)
is.numeric(summary_data_2$temperature)

lm<-lm(mean~temperature*flour,data=summary_data_2)
anova(lm)


