library(dplyr)
library(tidyr)
library(ggplot2)

# Simulation information
geneLength = 10000
postPolyALength	= 1800
promoterLength = 200
pauseSite = 40
resultFile = "result/default.result.txt"

# Read result file
res = read.table(resultFile,
	col.names = c("time", "start", "end", "count"))

# Bin 100 bases
res100 = res %>%
	mutate(pos = floor(start / 100)*100 - promoterLength) %>%
	group_by(time, pos) %>%
	summarize(count = sum(count))

# Plot
pdf("result/default.plot.pdf", width = 6, height = 4)
ggplot(data = res100, aes(x = pos, y = count)) +
	geom_bar(stat = "identity", aes(fill = time)) +
	facet_grid(time ~ .)
dev.off()

