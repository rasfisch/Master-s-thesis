library(tidyverse)

pc_k <- read_csv("/cloud/project/output_corrected/selection_model_sizes_corrected.csv") %>%
  select(date, pc_k)

summary(pc_k$pc_k)

pc_k %>%
  count(pc_k) %>%
  arrange(desc(n))

ggplot(pc_k, aes(x = date, y = pc_k)) +
  geom_line() +
  labs(
    title = "Number of principal components selected over time",
    x = NULL,
    y = "Selected number of PCs"
  ) +
  theme_minimal()