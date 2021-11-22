library(tidyverse)
library(readxl)
library(utils)
library(httr)
library(zoo)
library(cowplot)
library(purrr)
library(poweRlaw)
library(randomcoloR)
library(igraph)
library(Rcpp)
library(RcppArmadillo)
library(Matrix)
library(EpiEstim)
library(igraph)
library(treemapify)
library(gganimate)
library(rstan)
library(fixest)

#########################################################################################################
#                                              FUNCTIONS                                                #
#########################################################################################################

sourceCpp("RCPP/simulate_epidemic.cpp")

create_subpopulation <- function(adjacency_matrix, node_ids, degree_distribution) {
  nb_nodes <- length(node_ids)
  
  degrees <- degree_distribution(nb_nodes)
  
  # we need an even number of stubs to form edges
  while (sum(degrees) %% 2 != 0) {
    degrees <- degree_distribution(nb_nodes)
  }
  
  stubs <- unlist(imap(degrees, function(n, i) rep(node_ids[i], n)))
  nb_stubs <- length(stubs)
  
  sampled_stub_pairs <- sample(1:nb_stubs, nb_stubs)
  for (i in 1:(nb_stubs / 2)) {
    adjacency_matrix[stubs[sampled_stub_pairs[i * 2 - 1]], stubs[sampled_stub_pairs[i * 2]]] <- 1
    adjacency_matrix[stubs[sampled_stub_pairs[i * 2]], stubs[sampled_stub_pairs[i * 2 - 1]]] <- 1
  }
  
  adjacency_matrix
}

create_network <- function(
  nb_networks,
  nb_subpopulations_by_network_distribution,
  subpopulation_size_distribution,
  within_network_degree_distribution,
  between_network_degree_distribution,
  within_network_travel_prob,
  between_network_travel_prob
) {
  nb_subpopulations_by_network <- nb_subpopulations_by_network_distribution(nb_networks)
  
  nb_subpopulations <- sum(nb_subpopulations_by_network)
  
  subpopulations <- tibble(
    id = 1:nb_subpopulations, 
    network = unlist(map(1:nb_networks, function(i) rep(i, nb_subpopulations_by_network[i]))),
    population = subpopulation_size_distribution(nb_subpopulations)
  )
  
  adjacency_matrix <- Matrix(
    data = 0,
    nrow = nb_subpopulations,
    ncol = nb_subpopulations,
    sparse = TRUE
  )
  
  for (i in 1:nb_networks) {
    end_range <- sum(nb_subpopulations_by_network[1:i])
    start_range <- end_range - nb_subpopulations_by_network[i] + 1
    
    adjacency_matrix <- create_subpopulation(
      adjacency_matrix,
      start_range:end_range,
      within_network_degree_distribution
    )
  }
  
  # remove self-loops or sum(adjacency_matrix[i,]) will be less than 1 when the node i has a self-loop
  diag(adjacency_matrix) <- 0
  
  within_network_degrees <- rowSums(adjacency_matrix)
  for (i in 1:nb_subpopulations) {
    if (within_network_degrees[i] > 0) {
      adjacency_matrix[i,] <- adjacency_matrix[i,] * within_network_travel_prob
    }
  }
  
  between_network_degrees <- between_network_degree_distribution(nb_networks)
  
  # we need an even number of stubs to form edges
  while (sum(between_network_degrees) %% 2 != 0) {
    between_network_degrees <- between_network_degree_distribution(nb_networks)
  }
  
  stubs <- unlist(
    imap(
      between_network_degrees,
      function(n, i) {
        end_range <- sum(nb_subpopulations_by_network[1:i])
        start_range <- end_range - nb_subpopulations_by_network[i] + 1
        
        sample(
          start_range:end_range,
          n,
          replace = TRUE
        )
      }
    )
  )
  nb_stubs <- length(stubs)
  
  sampled_stub_pairs <- sample(1:nb_stubs, nb_stubs)
  for (i in 1:(nb_stubs / 2)) {
    adjacency_matrix[stubs[sampled_stub_pairs[i * 2 - 1]], stubs[sampled_stub_pairs[i * 2]]] <- 1
    adjacency_matrix[stubs[sampled_stub_pairs[i * 2]], stubs[sampled_stub_pairs[i * 2 - 1]]] <- 1
  }
  
  for (i in 1:nb_subpopulations) {
    adjacency_matrix[i,] <- ifelse(adjacency_matrix[i,] == 1, between_network_travel_prob, adjacency_matrix[i,])
  }
  
  # remove self-loops or sum(adjacency_matrix[i,]) will be less than 1 when the node i has a self-loop
  diag(adjacency_matrix) <- 0
  
  for (i in 1:nb_subpopulations) {
    travel_prob <- sum(adjacency_matrix[i,])
    adjacency_matrix[i, i] <- 1 - travel_prob
  }
  
  list(
    subpopulations = subpopulations,
    adjacency_matrix = adjacency_matrix
  )
}

dir.create("Figures")

# TODO: only seed the first variant for 100 days or so to reflect the fact that the pandemic slows down in most
# neighboring countries

#########################################################################################################
#                                    MODEL WITH HOMOGENEOUS MIXING                                      #
#########################################################################################################

set.seed(17)

simulation_length <- 130
end_seeding <- 30
seeding_rate <- 1
seeding_length <- 30
dispersion_factor <- 0.1
mean_gt <- 4.8
sd_gt <- 1.7

subpopulations_1 <- tibble(
  id = 1,
  population = 10e6
)

between_subpopulations_travel_probs_1 <- Matrix(
  data = 1,
  nrow = 1,
  ncol = 1,
  sparse = TRUE
)

seeding_from_abroad_by_destination_1 <- array(
  data = rpois(
    seeding_length,
    lambda = c(rep(seeding_rate, seeding_length), rep(0, simulation_length - seeding_length))
  ),
  dim = c(1, simulation_length, 1)
)

results_1 <- simulate_epidemic(
  subpopulations_1$population,
  simulation_length,
  array(
    data = 3,
    dim = c(1, 1, simulation_length)
  ),
  dispersion_factor,
  shape_gt = (mean_gt / sd_gt)^2,
  rate_gt = mean_gt / sd_gt^2,
  between_subpopulations_travel_probs_1,
  seeding_from_abroad_by_destination_1
)

results_1 <- tibble(
  t = seq_len(simulation_length),
  infections = results_1[1, , 1]
)

ggplot(results_1, aes(x = t, y = infections)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("Model of transmission with homogeneous mixing population") +
  xlab("Time") +
  ylab("Daily number of infections") +
  scale_y_continuous(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Model of transmission with homogeneous mixing population.png",
  width = 12,
  height = 6
)

set.seed(23)

simulation_length <- 250
start_lockdown <- 50
duration_lockdown <- 60

seeding_from_abroad_by_destination_2 <- array(
  data = rpois(
    seeding_length,
    lambda = c(rep(seeding_rate, seeding_length), rep(0, simulation_length - seeding_length))
  ),
  dim = c(1, simulation_length, 1)
)

results_2 <- simulate_epidemic(
  subpopulations_1$population,
  simulation_length,
  array(
    data = c(rep(3, start_lockdown - 1), rep(0.75, duration_lockdown), rep(3, simulation_length - duration_lockdown - start_lockdown + 1)),
    dim = c(1, 1, simulation_length)
  ),
  dispersion_factor,
  shape_gt = (mean_gt / sd_gt)^2,
  rate_gt = mean_gt / sd_gt^2,
  between_subpopulations_travel_probs_1,
  seeding_from_abroad_by_destination_2
)

results_2 <- tibble(
  t = seq_len(simulation_length),
  infections = results_2[1, , 1]
)

ggplot(results_2, aes(x = t, y = infections)) +
  geom_line(size = 1, color = "steelblue") +
  geom_vline(aes(xintercept = start_lockdown, linetype = "start_lockdown"), color = "darkred", size = 1) +
  geom_vline(aes(xintercept = start_lockdown + duration_lockdown, linetype = "end_lockdown"), color = "darkgreen", size = 1) +
  theme_bw() +
  ggtitle("Model of transmission with homogeneous mixing population and varying basic reproduction number") +
  xlab("Time") +
  ylab("Daily number of infections") +
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(
    name = "lines",
    values = c(
      "start_lockdown" = 2,
      "end_lockdown" = 2
    ),
    labels = c(
      "Start of lockdown",
      "End of lockdown"
    ),
    guide = guide_legend(
      title = "",
      override.aes = list(color = c("darkred", "darkgreen"))
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Model of transmission with homogeneous mixing population and varying basic reproduction number.png",
  width = 12,
  height = 6
)

#########################################################################################################
#                                        PLOTS OF REAL DATA                                             #
#########################################################################################################

google_mobility_data_url <- "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv"

google_mobility_data <- read_csv(google_mobility_data_url) %>%
  pivot_longer(
    cols = ends_with("baseline"),
    names_to = "type",
    values_to = "change"
  ) %>%
  mutate(
    type = case_when(
      type == "grocery_and_pharmacy_percent_change_from_baseline" ~ "Grocery and pharmacy",
      type == "parks_percent_change_from_baseline" ~ "Parks",
      type == "residential_percent_change_from_baseline" ~ "Residential",
      type == "retail_and_recreation_percent_change_from_baseline" ~ "Retail and recreation",
      type == "transit_stations_percent_change_from_baseline" ~ "Transit stations",
      type == "workplaces_percent_change_from_baseline" ~ "Workplaces"
    )
  ) %>%
  filter(type != "Parks") %>%
  select(date, country_region, sub_region_1, sub_region_2, change, type)

covid_data_sweden_url <- "https://www.arcgis.com/sharing/rest/content/items/b5e7488e117749c19881cce45db13f7e/data"

GET(
  covid_data_sweden_url,
  write_disk(tf <- tempfile(fileext = ".xlsx"))
)

data_icu_sweden <- read_excel(
  tf,
  sheet = "Antal intensivvårdade per dag"
) %>%
  rename(
    date = Datum_vårdstart,
    icu = Antal_intensivvårdade
  ) %>%
  mutate(
    date = as.Date(date),
    rolling_average_icu = rollmean(icu, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  filter(
    date >= as.Date("2020-03-06"),
    date <= as.Date("2020-05-31")
  ) %>%
  select(
    date,
    rolling_average_icu
  )

unlink(tf)

google_mobility_data_sweden <- google_mobility_data %>%
  filter(
    country_region == "Sweden",
    is.na(sub_region_1)
  ) %>%
  group_by(type) %>%
  mutate(
    smoothed_change_from_baseline = rollmean(change, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  ungroup() %>%
  filter(
    date >= as.Date("2020-03-06"),
    date <= as.Date("2020-05-31")
  ) %>%
  select(
    date,
    type,
    smoothed_change_from_baseline
  )

mobility_ICU_sweden_plot_title <- ggdraw() + 
  draw_label(
    "Mobility and daily number of ICU admissions in Sweden during the spring of 2020",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

google_mobility_sweden_plot <- ggplot(google_mobility_data_sweden, aes(x = date, y = smoothed_change_from_baseline, group = type, color = type)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Change from baseline in mobility") +
  scale_color_discrete(name = "Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  theme(
    legend.position = "bottom"
  ) +
  labs(caption = "Source: Google - Chart by Philippe Lemoine (@phl43)")

icu_sweden_plot <- ggplot(data_icu_sweden, aes(x = date, y = rolling_average_icu)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Daily number of ICU admissions for COVID-19") +
  scale_y_continuous() +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source: Folkhälsomyndigheten - Chart by Philippe Lemoine (@phl43)")

plot_grid(
  mobility_ICU_sweden_plot_title,
  google_mobility_sweden_plot,
  icu_sweden_plot,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Mobility and ICU admissions in Sweden in the spring of 2020.png",
  width = 12,
  height = 12
)

google_mobility_data_florida <- google_mobility_data %>%
  filter(
    country_region == "United States",
    sub_region_1 == "Florida",
    is.na(sub_region_2)
  ) %>%
  group_by(type) %>%
  mutate(
    smoothed_change_from_baseline = rollmean(change, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  ungroup() %>%
  filter(
    date >= as.Date("2020-06-01"),
    date <= as.Date("2020-08-31")
  ) %>%
  select(
    date,
    type,
    smoothed_change_from_baseline
  )

# https://github.com/nytimes/covid-19-data
us_data_url <- "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"

incidence_data_florida <- read_csv(us_data_url) %>%
  filter(state == "Florida") %>%
  mutate(
    new_cases = cases - lag(cases, default = 0),
    new_cases_smoothed = rollmean(new_cases, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  filter(
    date >= as.Date("2020-06-01"),
    date <= as.Date("2020-08-31")
  ) %>%
  select(
    date,
    new_cases_smoothed
  )

mobility_incidence_florida_plot_title <- ggdraw() + 
  draw_label(
    "Mobility and daily number of COVID-19 cases in Florida during the summer of 2020",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

google_mobility_plot_florida <- ggplot(google_mobility_data_florida, aes(x = date, y = smoothed_change_from_baseline, group = type, color = type)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Change from baseline") +
  scale_color_discrete(name = "Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  theme(
    legend.position = "bottom"
  ) +
  labs(caption = "Source: Google - Chart by Philippe Lemoine (@phl43)")

incidence_plot_florida <- ggplot(incidence_data_florida, aes(x = date, y = new_cases_smoothed)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Daily number of COVID-19 cases") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source: New York Times - Chart by Philippe Lemoine (@phl43)")

plot_grid(
  mobility_incidence_florida_plot_title,
  google_mobility_plot_florida,
  incidence_plot_florida,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Mobility and incidence in Florida in the summer of 2020.png",
  width = 12,
  height = 12
)

google_mobility_data_madrid <- google_mobility_data %>%
  filter(
    country_region == "Spain",
    sub_region_1 == "Community of Madrid"
  ) %>%
  group_by(type) %>%
  mutate(
    smoothed_change_from_baseline = rollmean(change, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  ungroup() %>%
  filter(
    date >= as.Date("2020-12-10"),
    date <= as.Date("2021-02-28")
  ) %>%
  select(
    date,
    type,
    smoothed_change_from_baseline
  )

# https://cnecovid.isciii.es/covid19/#documentación-y-datos
spain_incidence_data_url <- "https://cnecovid.isciii.es/covid19/resources/casos_tecnica_ccaa.csv"

incidence_data_madrid <- read_csv(spain_incidence_data_url) %>%
  filter(ccaa_iso == "MD") %>%
  rename(
    date = fecha,
    new_cases = num_casos
  ) %>%
  mutate(
    new_cases_smoothed = rollmean(new_cases, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  filter(
    date >= as.Date("2020-12-10"),
    date <= as.Date("2021-02-28")
  )

mobility_incidence_madrid_plot_title <- ggdraw() + 
  draw_label(
    "Mobility and daily number of COVID-19 cases in Madrid during the winter of 2020/2021",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

google_mobility_plot_madrid <- ggplot(google_mobility_data_madrid, aes(x = date, y = smoothed_change_from_baseline, group = type, color = type)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Change from baseline") +
  scale_color_discrete(name = "Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  theme(
    legend.position = "bottom"
  ) +
  labs(caption = "Source: Google - Chart by Philippe Lemoine (@phl43)")

incidence_plot_madrid <- ggplot(incidence_data_madrid, aes(x = date, y = new_cases_smoothed)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Daily number of COVID-19 cases") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source: Centro Nacional de Epidemiología - Chart by Philippe Lemoine (@phl43)")

plot_grid(
  mobility_incidence_madrid_plot_title,
  google_mobility_plot_madrid,
  incidence_plot_madrid,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Mobility and incidence in Madrid in the winter of 2020-2021.png",
  width = 12,
  height = 12
)

google_mobility_data_france <- google_mobility_data %>%
  filter(
    country_region == "France",
    is.na(sub_region_1)
  ) %>%
  group_by(type) %>%
  mutate(
    smoothed_change_from_baseline = rollmean(change, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  ungroup() %>%
  filter(
    date >= as.Date("2020-09-01"),
    date <= as.Date("2020-10-31")
  ) %>%
  select(
    date,
    type,
    smoothed_change_from_baseline
  )

# https://www.data.gouv.fr/en/datasets/donnees-relatives-aux-resultats-des-tests-virologiques-covid-19/
incidence_data_france_url <- "https://www.data.gouv.fr/fr/datasets/r/dd0de5d9-b5a5-4503-930a-7b08dc0adc7c"

incidence_data_france <- read_delim(incidence_data_france_url, delim = ";") %>%
  filter(cl_age90 == "0") %>%
  mutate(
    date = as.Date(jour)
  ) %>%
  rename(
    cases = P
  ) %>%
  mutate(
    new_cases_smoothed = round(rollmean(cases, 7, fill = c(0, 0, 0), align = "right"))
  ) %>%
  filter(
    date >= as.Date("2020-09-01"),
    date <= as.Date("2020-10-31")
  ) %>%
  select(
    date,
    new_cases_smoothed
  )

# https://www.data.gouv.fr/en/datasets/indicateurs-de-suivi-de-lepidemie-de-covid-19/
reproduction_number_data_france_url <- "https://www.data.gouv.fr/en/datasets/r/381a9472-ce83-407d-9a64-1b8c23af83df"

reproduction_number_data_france <- read_csv(reproduction_number_data_france_url) %>%
  rename(
    date = extract_date,
    reproduction_number = R
  ) %>%
  filter(
    date >= as.Date("2020-09-01"),
    date <= as.Date("2020-10-31")
  ) %>%
  select(
    date,
    reproduction_number
  )

mobility_incidence_france_plot_title <- ggdraw() + 
  draw_label(
    "Mobility and daily number of COVID-19 cases in France during the fall of 2020",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

mobility_R_france_plot_title <- ggdraw() + 
  draw_label(
    "Mobility and effective reproduction number in France during the fall of 2020",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

google_mobility_plot_france <- ggplot(google_mobility_data_france, aes(x = date, y = smoothed_change_from_baseline, group = type, color = type)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Change from baseline") +
  scale_color_discrete(name = "Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  theme(
    legend.position = "bottom"
  ) +
  labs(caption = "Source: Google - Chart by Philippe Lemoine (@phl43)")

incidence_plot_france <- ggplot(incidence_data_france, aes(x = date, y = new_cases_smoothed)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Daily number of COVID-19 cases") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source: Santé publique France - Chart by Philippe Lemoine (@phl43)")

reproduction_number_plot_france <- ggplot(reproduction_number_data_france, aes(x = date, y = reproduction_number)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Date") +
  ylab("Effective reproduction number") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source: Santé publique France - Chart by Philippe Lemoine (@phl43)")

plot_grid(
  mobility_incidence_france_plot_title,
  google_mobility_plot_france,
  incidence_plot_france,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Mobility and incidence in France in the fall of 2020.png",
  width = 12,
  height = 12
)

plot_grid(
  mobility_R_france_plot_title,
  google_mobility_plot_france,
  reproduction_number_plot_france,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Mobility and effective reproduction number in France in the fall of 2020.png",
  width = 12,
  height = 12
)

overseas_departments <- c(
  "971",
  "972",
  "973",
  "974",
  "975",
  "976",
  "977",
  "978"
)

# source: https://www.insee.fr/fr/statistiques/4265429?sommaire=4265511
french_departments <- read_delim(
  "Data/départements_français_métropole.csv",
  delim = ";"
) %>%
  filter(!(code %in% overseas_departments))

# https://www.data.gouv.fr/en/datasets/donnees-relatives-aux-resultats-des-tests-virologiques-covid-19/
incidence_data_by_department_france_url <- "https://www.data.gouv.fr/en/datasets/r/406c6a23-e283-4300-9484-54e78c8ae675"

incidence_data_by_department_france <- read_delim(incidence_data_by_department_france_url, delim = ";") %>%
  filter(cl_age90 == "0") %>%
  inner_join(french_departments, by = c("dep" = "code")) %>%
  mutate(
    date = as.Date(jour)
  ) %>%
  rename(
    cases = P
  ) %>%
  group_by(dep) %>%
  mutate(
    new_cases_smoothed = round(rollmean(cases, 7, fill = c(0, 0, 0), align = "right")),
    new_cases_smoothed_normalized = new_cases_smoothed / max(new_cases_smoothed)
  ) %>%
  ungroup() %>%
  select(
    date,
    dep,
    name,
    cases,
    new_cases_smoothed_normalized
  )

ggplot(incidence_data_by_department_france, aes(x = date, y = name)) +
  geom_tile(aes(fill = new_cases_smoothed_normalized)) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle(
    "Daily number of cases by date of specimen for departments in France (normalized by the maximum value during the period in each location)"
  ) +
  xlab("Date") +
  ylab("Department") +
  scale_x_date(
    labels = scales::date_format("%Y/%m/%d"),
    date_breaks = "1 month"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    fill = "Normalized\ndaily number\nof cases",
    caption = "Source: Santé publique France - Chart by Philippe Lemoine (@phl43)"
  )

ggsave(
  "Figures/Daily number of cases by date of specimen for departments in France.png",
  width = 24,
  height = 24
)

google_mobility_data_england <- google_mobility_data %>%
  filter(
    country_region == "United Kingdom",
    is.na(sub_region_1)
  ) %>%
  group_by(type) %>%
  mutate(
    smoothed_change_from_baseline = rollmean(change, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  ungroup() %>%
  filter(
    date >= as.Date("2021-05-15"),
    date <= as.Date("2021-11-10")
  ) %>%
  select(
    date,
    type,
    smoothed_change_from_baseline
  )

incidence_data_england_url <- "https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDate&format=csv"

incidence_data_england <- read_csv(incidence_data_england_url) %>%
  mutate(
    new_cases_smoothed = rollmean(newCasesBySpecimenDate, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  filter(
    date >= as.Date("2021-05-15"),
    date <= as.Date("2021-11-10")
  ) %>%
  select(
    date,
    new_cases_smoothed
  )

ggplot(incidence_data_england, aes(x = date, y = new_cases_smoothed)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("Daily number of COVID-19 cases in England by date of specimen (7-day moving average)") +
  xlab("Date") +
  ylab("Daily number of COVID-19 cases") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "15 day"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Public Health England - Chart by Philippe Lemoine (@phl43)")

ggsave(
  "Figures/Incidence in England in the summer and fall of 2021.png",
  width = 12,
  height = 6
)

ggplot(google_mobility_data_england, aes(x = date, y = smoothed_change_from_baseline, group = type, color = type)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Change from baseline of Google mobility data in England (7-day moving average)") +
  xlab("Date") +
  ylab("Change from baseline") +
  scale_color_discrete(name = "Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "15 day"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(caption = "Source: Google - Chart by Philippe Lemoine (@phl43)")

ggsave(
  "Figures/Mobility in England in the summer and fall of 2021.png",
  width = 12,
  height = 6
)

# source: https://www.data.gouv.fr/en/datasets/donnees-de-laboratoires-pour-le-depistage-indicateurs-sur-les-mutations/
url_department_mutations_data <- "https://www.data.gouv.fr/en/datasets/r/4d3e5a8b-9649-4c41-86ec-5420eb6b530c"

prevalence_mutations_by_department <- read_delim(url_department_mutations_data, delim = ";") %>%
  filter(!(dep %in% overseas_departments)) %>%
  mutate(
    # the prevalence of Delta in the dataset is the average for d to d + 7, so
    # I assign the estimate for a given week to d + 4
    date = ymd(str_extract(semaine, "^[:digit:]{4}-[:digit:]{2}-[:digit:]{2}")) + 3,
    nb_crib_E484K = nb_A0 + nb_A1,
    prevalence_E484K = ifelse(nb_crib_E484K > 0, nb_A1 / nb_crib_E484K, NA),
    prevalence_non_E484K = 1 - prevalence_E484K,
    nb_crib_E484Q = nb_B0 + nb_B1,
    prevalence_E484Q = ifelse(nb_crib_E484Q > 0, nb_B1 / nb_crib_E484Q, NA),
    prevalence_non_E484Q = 1 - prevalence_E484Q,
    nb_crib_L452R = nb_C0 + nb_C1,
    prevalence_L452R = ifelse(nb_crib_L452R > 0, nb_C1 / nb_crib_L452R, NA),
    prevalence_non_L452R = 1 - prevalence_L452R
  ) %>%
  select(
    dep,
    date,
    prevalence_E484K,
    prevalence_non_E484K,
    prevalence_E484Q,
    prevalence_non_E484Q,
    prevalence_L452R,
    prevalence_non_L452R,
    nb_crib_E484K,
    nb_crib_E484Q,
    nb_crib_L452R
  )

mean_gt <- 4.8

department_combined_data <- incidence_data_by_department_france %>%
  group_by(dep) %>%
  mutate(
    weekly_cases_total = rollsum(cases, 7, fill = rep(0, 6), align = "left"),
    weekly_growth_factor_total = weekly_cases_total / lag(weekly_cases_total, 7),
    R_total = weekly_growth_factor_total ^ (mean_gt / 7)
  ) %>%
  inner_join(prevalence_mutations_by_department %>% mutate(date = date - 3), by = c("dep", "date")) %>%
  mutate(
    weekly_cases_E484K = weekly_cases_total * prevalence_E484K,
    weekly_cases_E484Q = weekly_cases_total * prevalence_E484Q,
    weekly_cases_L452R = weekly_cases_total * prevalence_L452R,
    weekly_cases_non_E484K = weekly_cases_total * prevalence_non_E484K,
    weekly_cases_non_E484Q = weekly_cases_total * prevalence_non_E484Q,
    weekly_cases_non_L452R = weekly_cases_total * prevalence_non_L452R,
    weekly_growth_factor_E484K = weekly_cases_E484K / lag(weekly_cases_E484K, 7),
    weekly_growth_factor_E484Q = weekly_cases_E484Q / lag(weekly_cases_E484Q, 7),
    weekly_growth_factor_L452R = weekly_cases_L452R / lag(weekly_cases_L452R, 7),
    weekly_growth_factor_non_E484K = weekly_cases_non_E484K / lag(weekly_cases_non_E484K, 7),
    weekly_growth_factor_non_E484Q = weekly_cases_non_E484Q / lag(weekly_cases_non_E484Q, 7),
    weekly_growth_factor_non_L452R = weekly_cases_non_L452R / lag(weekly_cases_non_L452R, 7),
    R_E484K = weekly_growth_factor_E484K ^ (mean_gt / 7),
    R_E484Q = weekly_growth_factor_E484Q ^ (mean_gt / 7),
    R_L452R = weekly_growth_factor_L452R ^ (mean_gt / 7),
    R_non_E484K = weekly_growth_factor_non_E484K ^ (mean_gt / 7),
    R_non_E484Q = weekly_growth_factor_non_E484Q ^ (mean_gt / 7),
    R_non_L452R = weekly_growth_factor_non_L452R ^ (mean_gt / 7),
    advantage_E484K = R_E484K / R_non_E484K - 1,
    advantage_E484Q = R_E484Q / R_non_E484Q - 1,
    advantage_L452R = R_L452R / R_non_L452R - 1,
    prevalence_E484K = lag(prevalence_E484K, 7),
    prevalence_E484Q = lag(prevalence_E484Q, 7),
    prevalence_L452R = lag(prevalence_L452R, 7)
  ) %>%
  ungroup() %>%
  filter(
    date %in% seq(ymd("2021-05-31"), ymd("2021-09-27"), by = "7 day") &
      nb_crib_L452R >= 30
  ) %>%
  mutate(
    period = paste0("Week ", as.integer(date - ymd("2021-01-04")) / 7, " to week ", as.integer(date - ymd("2021-01-04")) / 7 + 1),
    t = as.integer(date - ymd("2021-01-04")) / 7 - 21
  ) %>%
  na.omit() %>%
  arrange(period) %>%
  select(
    t,
    period,
    dep,
    weekly_cases_L452R,
    weekly_cases_non_L452R,
    weekly_growth_factor_L452R,
    weekly_growth_factor_non_L452R,
    R_total,
    R_L452R,
    R_non_L452R,
    advantage_L452R,
    prevalence_L452R
  )

department_combined_data$period <- fct_relevel(
  department_combined_data$period,
  c(
    "Week 22 to week 23",
    "Week 23 to week 24",
    "Week 24 to week 25",
    "Week 25 to week 26",
    "Week 26 to week 27",
    "Week 27 to week 28",
    "Week 28 to week 29",
    "Week 29 to week 30",
    "Week 30 to week 31",
    "Week 31 to week 32",
    "Week 32 to week 33",
    "Week 33 to week 34",
    "Week 34 to week 35",
    "Week 35 to week 36",
    "Week 36 to week 37",
    "Week 37 to week 38",
    "Week 38 to week 39"
  )
)

ggplot(department_combined_data, aes(x = prevalence_L452R, y = advantage_L452R)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.5, 3)) +
  theme_bw() +
  ggtitle(
    "Delta's transmission advantage vs. prevalence of Delta at the department level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Delta's transmission advantage") +
  scale_color_discrete(name = "Period") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave(
  "Figures/Delta's transmission advantage vs. prevalence of Delta at the department level in France.png",
  width = 12,
  height = 6
)

#########################################################################################################
#                                     MODEL WITH POPULATION STRUCTURE                                   #
#########################################################################################################

set.seed(71)

nb_networks <- 100
mean_nb_subpopulations_by_network <- 100
min_subpopulation_population_size <- 250
alpha_subpopulation_population_size <- 2.9
mean_between_network_degree <- 5
mean_within_network_degree_distribution <- 5
simulation_length <- 365
within_network_travel_prob <- 0.05
between_network_travel_prob <- 0.0001
mean_R0 <- 2.5
sd_R0 <- 0.1
seeding_rate <- 1
dispersion_factor <- 0.1
mean_gt <- 4.8
sd_gt <- 1.7

nb_subpopulations_by_network_distribution <- partial(
  rpois,
  lambda = mean_nb_subpopulations_by_network
)

subpopulation_size_distribution <- partial(
  rpldis,
  xmin = min_subpopulation_population_size,
  alpha = alpha_subpopulation_population_size
)

subpopulation_R0_generation_function <- function(nb_subpopulations, nb_variants, simulation_length) {
  R0_by_subpopulation <- rnorm(nb_subpopulations, mean = mean_R0, sd = sd_R0)
  
  array(
    data = rep(R0_by_subpopulation, each = nb_variants * simulation_length),
    dim = c(nb_subpopulations, nb_variants, simulation_length)
  )
}

between_network_degree_distribution <- partial(rpois, lambda = mean_between_network_degree)
within_network_degree_distribution <- partial(rpois, lambda = mean_within_network_degree_distribution)

full_network_3 <- create_network(
  nb_networks,
  nb_subpopulations_by_network_distribution,
  subpopulation_size_distribution,
  within_network_degree_distribution,
  between_network_degree_distribution,
  within_network_travel_prob,
  between_network_travel_prob
)

subpopulations_3 <- full_network_3$subpopulations
between_subpopulations_travel_probs_3 <- full_network_3$adjacency_matrix

R0_3 <- subpopulation_R0_generation_function(nrow(subpopulations_3), 1, simulation_length)

seeding_from_abroad <- rpois(
  simulation_length,
  lambda = seeding_rate
)

seeding_from_abroad_by_destination_3 <- array(
  data = rep(0, nrow(subpopulations_3) * simulation_length),
  dim = c(nrow(subpopulations_3), simulation_length, 1)
)

# for some reason, it seems that sampling doesn't work properly in Rcpp, so I'm doing this part here since it's very fast
# in R anyway (see https://stackoverflow.com/questions/60119621/get-the-same-sample-of-integers-from-rcpp-as-base-r for
# a possible explanation, but I'm not sure what is discussed in the comments of that question can be what's going on here)
for (i in 1:simulation_length) {
  if (seeding_from_abroad[i] > 0) {
    destination_units <- sample(1:nrow(subpopulations_3), seeding_from_abroad[i], replace = TRUE)
    
    for (j in 1:seeding_from_abroad[i]) {
      seeding_from_abroad_by_destination_3[destination_units[j], i, 1] <- seeding_from_abroad_by_destination_3[destination_units[j], i, 1] + 1
    }
  }
}

results_3 <- simulate_epidemic(
  subpopulations_3$population,
  simulation_length,
  R0_3,
  dispersion_factor,
  shape_gt = (mean_gt / sd_gt)^2,
  rate_gt = mean_gt / sd_gt^2,
  between_subpopulations_travel_probs_3,
  seeding_from_abroad_by_destination_3
)

results_3 <- map2_dfr(
  rep(subpopulations_3$id, 1),
  rep(1L, nrow(subpopulations_3)),
  function(id, v) {
    tibble(
      t = seq_len(simulation_length),
      subpopulation = id,
      infections = results_3[id, , v],
      variant = v
    )
  }
) %>%
  inner_join(subpopulations_3, by = c("subpopulation" = "id"))

results_national_3 <- results_3 %>%
  group_by(t) %>%
  summarize(infections = sum(infections))

attack_rate_over_time_3 <- results_3 %>%
  group_by(subpopulation, population, t) %>%
  summarize(
    infections = sum(infections)
  ) %>%
  mutate(
    cumulative_infections = cumsum(infections)
  ) %>%
  ungroup() %>%
  group_by(t) %>%
  summarize(
    attack_rate = sum(cumulative_infections) / sum(population)
  )

start_estimation <- filter(results_national_3, infections > 100)$t[1]
end_estimation <- filter(results_national_3, t > start_estimation & infections < 100)$t[1]
end_estimation <- ifelse(is.na(end_estimation), simulation_length, end_estimation)

effective_reproduction_number_3 <- estimate_R(
  results_national_3$infections,
  method = "parametric_si",
  config = make_config(
    t_start = start_estimation:end_estimation,
    t_end = start_estimation:end_estimation,
    method = "parametric_si",
    si_parametric_distr = "G",
    mean_si = mean_gt,
    std_si = sd_gt
  )
)$R %>%
  mutate(
    t = t_start,
    R = `Mean(R)`
  )

mobility_incidence_simulation_plot_title_3 <- ggdraw() + 
  draw_label(
    "Mobility and daily number of infections for the population as a whole in a simulation with complex population structure",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

effective_reproduction_number_plot_3 <- ggplot(effective_reproduction_number_3, aes(x = t, y = R)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Time") +
  ylab("Effective reproduction number") +
  xlim(0, 365) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

national_incidence_plot_3 <- ggplot(results_national_3, aes(x = t, y = infections)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Time") +
  ylab("Daily number of infections") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

plot_grid(
  mobility_incidence_simulation_plot_title_3,
  effective_reproduction_number_plot_3,
  national_incidence_plot_3,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Effective reproduction number and daily number of infections for the population as a whole in a simulation with complex population structure.png",
  width = 12,
  height = 12
)

igraph_network <- graph_from_adjacency_matrix(
  between_subpopulations_travel_probs_3,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

color_palette <- randomColor(nb_networks, luminosity = "light")
color_by_subpopulation <- map_dfr(subpopulations_3$network, function(network) as_tibble(t(col2rgb(color_palette[network])))) %>%
  rename(r = red, g = green, b = blue)

gexf_network <- igraph.to.gexf(
  igraph_network,
  nodesAtt = subpopulations_3 %>% select(population),
  nodesVizAtt = tibble(
    color = color_by_subpopulation
  ),
  edgesVizAtt = tibble(
    color = tibble(r = 192, g = 192, b = 192)
  )
)
gexf_network$edgesVizAtt$color <- tibble(r = 192, g = 192, b = 192)
write.gexf(gexf_network, output = "network.gexf")

data_for_animation <- results_3 %>%
  group_by(subpopulation) %>%
  mutate(
    incidence_rate = infections / population * 1e5,
    cumulative_infections = cumsum(infections),
    attack_rate = cumulative_infections / population
  ) %>%
  ungroup()

# attack_rates_by_subpopulation_plot <- ggplot(data_for_animation, aes(area = population, fill = attack_rate, subgroup = network)) +
#   geom_treemap(layout = "fixed") +
#   transition_time(t) +
#   theme_bw() +
#   scale_fill_viridis_c(
#     name = "Attack rate",
#     labels = function(x) paste0(x * 100, "%"),
#     values = scales::rescale(c(0, 0.25, 0.75, 1)),
#     option = "inferno"
#   ) +
#   ggtitle(
#     "Attack rate in each subpopulation over the course of the simulation",
#     subtitle = "(t = {frame_time})"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5)
#   )
# 
# attack_rates_by_subpopulation_animation <- animate(
#   attack_rates_by_subpopulation_plot,
#   width = 1200,
#   height = 1200,
#   duration = 120,
#   fps = 3,
#   renderer = av_renderer()
# )
# 
# anim_save(
#   "How the virus spreads in a model with complex population structure.mp4",
#   attack_rates_by_subpopulation_animation
# )

nb_regions <- 15
network_geographical_concentration <- 10

nb_subpopulations_3 <- nrow(subpopulations_3)
nb_subpopulations_by_network_3 <- unlist(map(seq_len(nb_networks), function (network) sum(subpopulations_3$network == network)))
subpopulations_3$region <- rep(0L, nb_subpopulations_3)
main_region_by_network_3 <- sample(1:nb_regions, nb_networks, replace = TRUE)
for (i in 1:nb_networks) {
  end_range <- sum(nb_subpopulations_by_network_3[1:i])
  start_range <- end_range - nb_subpopulations_by_network_3[i] + 1
  subpopulations_3$region[start_range:end_range] <- sample(
    1:nb_regions,
    end_range - start_range + 1,
    prob = c(rep(1, main_region_by_network_3[i] - 1), network_geographical_concentration, rep(1, nb_regions - main_region_by_network_3[i])),
    replace = TRUE
  )
}

results_3 <- inner_join(
  results_3,
  subpopulations_3 %>% select(id, region),
  by = c("subpopulation" = "id")
)

results_regional_3 <- results_3 %>%
  group_by(t, region) %>%
  summarize(infections = sum(infections)) %>%
  ungroup()

normalized_incidence_by_region_3 <- results_regional_3 %>%
  group_by(region) %>%
  mutate(
    incidence_normalized = infections / max(infections)
  ) %>%
  ungroup()

ggplot(normalized_incidence_by_region_3, aes(x = t, y = factor(region))) +
  geom_tile(aes(fill = incidence_normalized), interpolate = TRUE) +
  scale_fill_viridis_c(
    #direction = -1
  ) +
  theme_bw() +
  ggtitle(
    "Daily number of infections by region in a simulation with complex population structure (normalized by the maximum value during the period in each location)"
  ) +
  xlab("Time") +
  ylab("Region") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    fill = "Normalized\ndaily number\nof cases"
  )

ggsave(
  "Figures/Daily number of infections by region in a simulation with complex population structure (normalized by the maximum value during the period in each location).png",
  width = 16,
  height = 12
)

#########################################################################################################
#          MODEL WITH POPULATION STRUCTURE BUT MORE CONNECTIVITY BETWEEN SUBNETWORKS                    #
#########################################################################################################

set.seed(37)

mean_between_network_degree <- 50
between_network_travel_prob <- 0.01

full_network_4 <- create_network(
  nb_networks,
  nb_subpopulations_by_network_distribution,
  subpopulation_size_distribution,
  within_network_degree_distribution,
  between_network_degree_distribution,
  within_network_travel_prob,
  between_network_travel_prob
)

subpopulations_4 <- full_network_4$subpopulations
between_subpopulations_travel_probs_4 <- full_network_4$adjacency_matrix

R0_4 <- subpopulation_R0_generation_function(nrow(subpopulations_4), 1, simulation_length)

seeding_from_abroad <- rpois(
  simulation_length,
  lambda = seeding_rate
)

seeding_from_abroad_by_destination_4 <- array(
  data = rep(0, nrow(subpopulations_4) * simulation_length),
  dim = c(nrow(subpopulations_4), simulation_length, 1)
)

# for some reason, it seems that sampling doesn't work properly in Rcpp, so I'm doing this part here since it's very fast
# in R anyway (see https://stackoverflow.com/questions/60119621/get-the-same-sample-of-integers-from-rcpp-as-base-r for
# a possible explanation, but I'm not sure what is discussed in the comments of that question can be what's going on here)
for (i in 1:simulation_length) {
  if (seeding_from_abroad[i] > 0) {
    destination_units <- sample(1:nrow(subpopulations_4), seeding_from_abroad[i], replace = TRUE)
    
    for (j in 1:seeding_from_abroad[i]) {
      seeding_from_abroad_by_destination_4[destination_units[j], i, 1] <- seeding_from_abroad_by_destination_4[destination_units[j], i, 1] + 1
    }
  }
}

results_4 <- simulate_epidemic(
  subpopulations_4$population,
  simulation_length,
  R0_4,
  dispersion_factor,
  shape_gt = (mean_gt / sd_gt)^2,
  rate_gt = mean_gt / sd_gt^2,
  between_subpopulations_travel_probs_4,
  seeding_from_abroad_by_destination_4
)

results_4 <- map2_dfr(
  rep(subpopulations_4$id, 1),
  rep(1L, nrow(subpopulations_4)),
  function(id, v) {
    tibble(
      t = seq_len(simulation_length),
      subpopulation = id,
      infections = results_4[id, , v],
      variant = v
    )
  }
) %>%
  inner_join(subpopulations_4, by = c("subpopulation" = "id"))

results_national_4 <- results_4 %>%
  group_by(t) %>%
  summarize(infections = sum(infections))

ggplot(results_national_4, aes(x = t, y = infections)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("Daily number of infections for the population as a whole in a simulation with complex population structure") +
  xlab("Time") +
  ylab("Daily number of infections") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Daily number of infections for the population as a whole in a simulation with complex population structure but more connectivity between subnetworks.png",
  width = 12,
  height = 6
)

attack_rate_over_time_4 <- results_4 %>%
  group_by(subpopulation, population, t) %>%
  summarize(
    infections = sum(infections)
  ) %>%
  mutate(
    cumulative_infections = cumsum(infections)
  ) %>%
  ungroup() %>%
  group_by(t) %>%
  summarize(
    attack_rate = sum(cumulative_infections) / sum(population)
  )

#########################################################################################################
#               MODEL WITH POPULATION STRUCTURE AND TWO EQUALLY TRANSMISSIBLE VARIANTS                  #
#########################################################################################################

set.seed(112)

mean_between_network_degree <- 5
between_network_travel_prob <- 0.0001

end_seeding_first_variant <- 150
start_seeding_second_variant <- 180

full_network_5 <- create_network(
  nb_networks,
  nb_subpopulations_by_network_distribution,
  subpopulation_size_distribution,
  within_network_degree_distribution,
  between_network_degree_distribution,
  within_network_travel_prob,
  between_network_travel_prob
)

subpopulations_5 <- full_network_5$subpopulations
between_subpopulations_travel_probs_5 <- full_network_5$adjacency_matrix

R0_5 <- subpopulation_R0_generation_function(nrow(subpopulations_5), 2, simulation_length)

seeding_from_abroad_variant1 <- rpois(
  simulation_length,
  lambda = c(rep(seeding_rate, end_seeding_first_variant), rep(0, simulation_length - end_seeding_first_variant))
)

seeding_from_abroad_variant2 <- rpois(
  simulation_length,
  lambda = c(rep(0, start_seeding_second_variant - 1), rep(seeding_rate, simulation_length - start_seeding_second_variant))
)

seeding_from_abroad_by_destination_5 <- array(
  data = rep(0, nrow(subpopulations_5) * simulation_length),
  dim = c(nrow(subpopulations_5), simulation_length, 2)
)

# for some reason, it seems that sampling doesn't work properly in Rcpp, so I'm doing this part here since it's very fast
# in R anyway (see https://stackoverflow.com/questions/60119621/get-the-same-sample-of-integers-from-rcpp-as-base-r for
# a possible explanation, but I'm not sure what is discussed in the comments of that question can be what's going on here)
for (i in 1:simulation_length) {
  if (seeding_from_abroad_variant1[i] > 0) {
    destination_units <- sample(1:nrow(subpopulations_5), seeding_from_abroad_variant1[i], replace = TRUE)
    
    for (j in 1:seeding_from_abroad_variant1[i]) {
      seeding_from_abroad_by_destination_5[destination_units[j], i, 1] <- seeding_from_abroad_by_destination_5[destination_units[j], i, 1] + 1
    }
  }
  
  if (seeding_from_abroad_variant2[i] > 0) {
    destination_units <- sample(1:nrow(subpopulations_5), seeding_from_abroad_variant2[i], replace = TRUE)
    
    for (j in 1:seeding_from_abroad_variant1[i]) {
      seeding_from_abroad_by_destination_5[destination_units[j], i, 2] <- seeding_from_abroad_by_destination_5[destination_units[j], i, 2] + 1
    }
  }
}

results_5 <- simulate_epidemic(
  subpopulations_5$population,
  simulation_length,
  R0_5,
  dispersion_factor,
  shape_gt = (mean_gt / sd_gt)^2,
  rate_gt = mean_gt / sd_gt^2,
  between_subpopulations_travel_probs_5,
  seeding_from_abroad_by_destination_5
)

results_5 <- map2_dfr(
  rep(subpopulations_5$id, 2),
  c(rep(1L, nrow(subpopulations_5)), rep(2L, nrow(subpopulations_5))),
  function(id, v) {
    tibble(
      t = seq_len(simulation_length),
      subpopulation = id,
      infections = results_5[id, , v],
      variant = v
    )
  }
) %>%
  inner_join(subpopulations_5, by = c("subpopulation" = "id"))

results_5$variant <- factor(results_5$variant)

results_national_5 <- results_5 %>%
  group_by(t) %>%
  summarize(infections = sum(infections))

results_national_by_variant <- results_5 %>%
  group_by(t, variant) %>%
  summarize(infections = sum(infections))

mobility_incidence_simulation_plot_title_5 <- ggdraw() + 
  draw_label(
    "Daily number of infections for the population as a whole in a simulation with\ncomplex population structure and two equally transmissible variants",
    fontface = "bold"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white")
  )

national_incidence_plot_5 <- ggplot(results_national_5, aes(x = t, y = infections)) +
  geom_line(size = 1, color = "steelblue") +
  theme_bw() +
  ggtitle("") +
  xlab("Time") +
  ylab("Daily number of infections") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

national_incidence_by_variant_plot <- ggplot(results_national_by_variant, aes(x = t, y = infections, group = variant, color = variant)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("") +
  xlab("Time") +
  ylab("Daily number of infections") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

plot_grid(
  mobility_incidence_simulation_plot_title_5,
  national_incidence_plot_5,
  national_incidence_by_variant_plot,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)

ggsave(
  "Figures/Daily number of infections for the population as a whole in a simulation with complex population structure and two equally transmissible variants.png",
  width = 12,
  height = 12
)

attack_rate_over_time_5 <- results_5 %>%
  group_by(subpopulation, population, t) %>%
  summarize(
    infections = sum(infections)
  ) %>%
  mutate(
    cumulative_infections = cumsum(infections)
  ) %>%
  ungroup() %>%
  group_by(t) %>%
  summarize(
    attack_rate = sum(cumulative_infections) / sum(population)
  )

results_national_by_variant <- results_national_by_variant %>%
  pivot_wider(
    names_from = variant,
    names_prefix = "infections_variant",
    values_from = infections
  )

effective_reproduction_number_variant1 <- estimate_R(
  results_national_by_variant$infections_variant1,
  method = "parametric_si",
  config = make_config(
    t_start = 2:(simulation_length - 6),
    t_end = (2 + 6):simulation_length,
    method = "parametric_si",
    si_parametric_distr = "G",
    mean_si = mean_gt,
    std_si = sd_gt
  )
)$R %>%
  mutate(
    variant = 1L,
    t = t_start,
    R = `Mean(R)`
  ) %>%
  select(variant, t, R)

effective_reproduction_number_variant2 <- estimate_R(
  results_national_by_variant$infections_variant2,
  method = "parametric_si",
  config = make_config(
    t_start = 2:(simulation_length - 6),
    t_end = (2 + 6):simulation_length,
    method = "parametric_si",
    si_parametric_distr = "G",
    mean_si = mean_gt,
    std_si = sd_gt
  )
)$R %>%
  mutate(
    variant = 2L,
    t = t_start,
    R = `Mean(R)`
  ) %>%
  select(variant, t, R)

effective_reproduction_number_5 <- bind_rows(effective_reproduction_number_variant1, effective_reproduction_number_variant2) %>%
  inner_join(results_national_by_variant, by = "t") %>%
  filter(infections_variant1 >= 100 & infections_variant2 >= 100)

transmission_advantage <- effective_reproduction_number_5 %>%
  pivot_wider(
    names_from = variant,
    names_prefix = "R_variant",
    values_from = R
  ) %>%
  mutate(
    advantage = R_variant2 / R_variant1 - 1
  )

ggplot(transmission_advantage, aes(x = t, y = advantage)) +
  geom_line(size = 1, color = "steelblue") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_bw() +
  ggtitle("Second variant's transmission advantage in a simulation with complex population structure and two equally transmissible variants") +
  xlab("Time") +
  ylab("Transmission advantage") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Second variant's transmission advantage in a simulation with complex population structure and two equally transmissible variants.png",
  width = 12,
  height = 6
)

nb_subpopulations <- nrow(subpopulations_5)
nb_subpopulations_by_network <- unlist(map(seq_len(nb_networks), function (network) sum(subpopulations_5$network == network)))
subpopulations_5$region <- rep(0L, nb_subpopulations)
main_region_by_network <- sample(1:nb_regions, nb_networks, replace = TRUE)
for (i in 1:nb_networks) {
  end_range <- sum(nb_subpopulations_by_network[1:i])
  start_range <- end_range - nb_subpopulations_by_network[i] + 1
  subpopulations_5$region[start_range:end_range] <- sample(
    1:nb_regions,
    end_range - start_range + 1,
    prob = c(rep(1, main_region_by_network[i] - 1), network_geographical_concentration, rep(1, nb_regions - main_region_by_network[i])),
    replace = TRUE
  )
}

results_5 <- inner_join(
  results_5,
  subpopulations_5 %>% select(id, region),
  by = c("subpopulation" = "id")
)

results_by_region_and_variant <- results_5 %>%
  group_by(t, region, variant) %>%
  summarize(
    infections = sum(infections)
  ) %>%
  pivot_wider(
    names_from = variant,
    names_prefix = "infections_variant",
    values_from = infections
  )

effective_reproduction_number_by_region <- tibble()

for (r in seq_len(nb_regions)) {
  results_by_variant <- results_by_region_and_variant %>%
    filter(region == r)
  
  effective_reproduction_number_variant1 <- estimate_R(
    results_by_variant$infections_variant1,
    method = "parametric_si",
    config = make_config(
      t_start = 2:(simulation_length - 6),
      t_end = (2 + 6):simulation_length,
      method = "parametric_si",
      si_parametric_distr = "G",
      mean_si = mean_gt,
      std_si = sd_gt
    )
  )$R %>%
    mutate(
      region = r,
      variant = 1L,
      t = t_start,
      R = `Mean(R)`,
    ) %>%
    select(region, variant, t, R)
  
  effective_reproduction_number_variant2 <- estimate_R(
    results_by_variant$infections_variant2,
    method = "parametric_si",
    config = make_config(
      t_start = 2:(simulation_length - 6),
      t_end = (2 + 6):simulation_length,
      method = "parametric_si",
      si_parametric_distr = "G",
      mean_si = mean_gt,
      std_si = sd_gt
    )
  )$R %>%
    mutate(
      region = r,
      variant = 2L,
      t = t_start,
      R = `Mean(R)`,
    ) %>%
    select(region, variant, t, R)
  
  effective_reproduction_number_by_region <- bind_rows(effective_reproduction_number_variant1, effective_reproduction_number_variant2) %>%
    inner_join(results_by_region_and_variant, by = c("t", "region")) %>%
    filter(infections_variant1 >= 25 & infections_variant2 >= 25) %>%
    bind_rows(effective_reproduction_number_by_region) %>%
    select(region, t, variant, R)
}

transmission_advantage_by_region <- effective_reproduction_number_by_region %>%
  pivot_wider(
    names_from = variant,
    names_prefix = "R_variant",
    values_from = R
  ) %>%
  mutate(
    advantage = R_variant2 / R_variant1 - 1
  ) %>%
  filter(
    region %in% (effective_reproduction_number_by_region %>% group_by(region) %>% summarize(n = n()) %>% filter(n >= 20))$region
  )

ggplot(transmission_advantage_by_region, aes(x = t, y = advantage)) +
  geom_line(size = 1, color = "steelblue") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_bw() +
  ggtitle("Second variant's transmission advantage in a simulation with complex population structure and two equally transmissible variants") +
  xlab("Time") +
  ylab("Transmission advantage") +
  facet_wrap(~ region) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Second variant's transmission advantage by region in a simulation with complex population structure and two equally transmissible variants.png",
  width = 12,
  height = 6
)

#########################################################################################################
#         USING A STANDARD EPIDEMIOLOGICAL MODEL TO INFER THE EFFECT OF A FICTITIOUS NPI                #
#########################################################################################################

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

model <- stan_model("Stan/model.stan")

data_length <- 200
seeding_length <- 30
start_lockdown <- 105
end_lockdown <- 180

ggplot(results_national_3 %>% filter(t <= data_length), aes(x = t, y = infections)) +
  geom_line(size = 1, color = "steelblue") +
  annotate("rect", xmin = start_lockdown, xmax = end_lockdown, ymin = 0, ymax = Inf, alpha = 0.2, fill = "red") +
  theme_bw() +
  ggtitle(
    "Daily number of infections for the population as a whole in a simulation with complex population structure",
    subtitle = "(the model was told that a lockdown was in effect during the period indicated as shaded area)"
  ) +
  xlab("Time") +
  ylab("Daily number of infections") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Simulated data fed to a standard epidemiological model with homogeneous population mixing.png",
  width = 12,
  height = 6
)

X <- rollmean(
  c(rep(0, start_lockdown - 1), rep(1, end_lockdown - start_lockdown), rep(0, data_length - end_lockdown + 1)),
  7,
  fill = rep("extend", 7 - 1),
  align = "right"
)

stan_data <- list(
  N = data_length,
  population = sum(subpopulations_3$population),
  seeding_length = 30,
  prediction_length = data_length - seeding_length,
  infections = results_national_3$infections[1:data_length],
  X = X,
  mu = mean_gt,
  sigma = sd_gt
)

fit_mcmc <- sampling(
  model,
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 17,
  control = list(adapt_delta = 0.95, max_treedepth = 10)
)

output <- rstan::extract(fit_mcmc)
estimate_lockdown_effect <- 1 - exp(-mean(output$gamma))
quantile(1 - exp(-output$gamma), probs = c(0.1, 0.9))

#########################################################################################################
#               USING A SIMPLE ECONOMETRIC MODEL TO INFER THE EFFECT OF A FICTITIOUS NPI                #
#########################################################################################################

set.seed(51)

n_simulations <- 1000

lockdown_starting_range <- 100:140
lockdown_length <- 60

difflag <- 1

model_fits <- list()
for (i in 1:n_simulations) {
  lockdown_starting_times <- sample(lockdown_starting_range, nb_regions, replace = TRUE)
  
  data_for_regression <- results_regional_3 %>%
    filter(t <= data_length) %>%
    group_by(region) %>%
    group_modify(
      ~ tibble(
        t = .x$t,
        dc = rollsum(.x$infections, difflag, fill = rep(0, difflag - 1), align = "right"),
        logdc = ifelse(dc > 0, log(dc), -1),
        dlogdc = logdc - dplyr::lag(logdc, difflag, default = -1),
        lockdown = rollmean(
          c(rep(0, lockdown_starting_times[.y$region] - 1), rep(1, lockdown_length), rep(0, data_length - lockdown_starting_times[.y$region] - lockdown_length + 1)),
          difflag,
          fill = rep("extend", difflag - 1),
          align = "right"
        )
      )
    ) %>%
    ungroup()
  
  model_fits[[i]] <- feols(dlogdc ~ l(lockdown, 5) | region, data_for_regression, panel.id = ~region+t)
}

estimates <- map_dfr(
  1:n_simulations,
  function(i) tibble(
    coefficient = coef(model_fits[[i]]),
    se = se(model_fits[[i]]),
    lower = coefficient - 1.96 * se,
    upper = coefficient + 1.96 * se
  )
)

ggplot(estimates, aes(x = coefficient)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "steelblue") +
  theme_bw() +
  ggtitle("Estimates of the effect of a non-existent lockdown on the daily growth of infections in data simulated\nby a model assuming complex population structure according to a simple econometric model") +
  xlab("Estimate") +
  ylab("Density") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "Figures/Estimates of the effect of a non-existent lockdown on the daily growth of infections in data simulated with a model assuming complex population structure\naccording to a simple econometric model.png",
  width = 12,
  height = 6
)
