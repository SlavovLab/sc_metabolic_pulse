library(dplyr)
library(ggplot2)


#####################
# Figure 1 Spectra plot
#####################


plot_spectra <- read.csv('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/Users/aleduc/LYSF.csv')
plot_spectra$mz_values <- round(plot_spectra$mz_values,2)
plot_spectra <- plot_spectra %>% group_by(mz_values) %>% summarise(intensity_values = sum(intensity_values))

ggplot(plot_spectra,aes(x = mz_values,y = intensity_values)) + geom_line()

'LYSEFLGK'

spects_want <- c(618.85,
                 618.85 + .5,
                 618.85 + 1,
                 621.86,
                 621.86 +.5,
                 622.79,
                 626.86,
                 626.86+.5,
                 626.86+1,
                 629.87,
                 629.87+.5,
                 629.77+1)

plot_spectra_fuzz <- plot_spectra %>% filter(intensity_values < 5809)
plot_spectra_fuzz$type <- 'Fuzz'
plot_spectra_fuzz <- plot_spectra_fuzz %>% filter(!mz_values %in% spects_want ) 
#plot_spectra_fuzz <- plot_spectra_fuzz[sample(1:nrow(plot_spectra_fuzz),nrow(plot_spectra_fuzz)/10),]

plot_spectra2 <- plot_spectra %>% filter(mz_values %in% spects_want)
plot_spectra2$type <- 'pep'

plot_spectra2$type[plot_spectra2$mz_values < 620] <- '1'
plot_spectra2$type[plot_spectra2$mz_values > 620 & plot_spectra2$mz_values < 625] <- '2'
plot_spectra2$type[plot_spectra2$mz_values > 625 & plot_spectra2$mz_values < 629] <- '3'
plot_spectra2$type[plot_spectra2$mz_values > 629] <- '4'

plot_spectra3 <- rbind(plot_spectra2,plot_spectra_fuzz)

ggplot(plot_spectra3,aes(x = mz_values,y = intensity_values,fill = type)) + 
  geom_bar(stat = 'identity',width=.15) + theme_classic() + xlim(c(615,635)) +
  scale_fill_manual(values = c('#F6AA4C','#D27628','#F6AA4C','#D27628','gray85'))


ggplot() +
  ## wide bars for A–D
  geom_bar(
    data = subset(plot_spectra3, type %in% c("1", "2", "3", "4")),
    aes(x = mz_values, y = intensity_values, fill = type),
    stat  = "identity",
    width = 0.15
  ) +
  ## narrow bars for E
  geom_bar(
    data = subset(plot_spectra3, type == "Fuzz"),
    aes(x = mz_values, y = intensity_values, fill = type),
    stat  = "identity",
    width = 0.05
  ) +
  theme_classic() +
  xlim(615, 635) +
  scale_fill_manual(
    values = c(
      "1" = "#F6AA4C",
      "2" = "#D27628",
      "3" = "#F6AA4C",
      "4" = "#D27628",
      "Fuzz" = "gray85"
    )
  )


