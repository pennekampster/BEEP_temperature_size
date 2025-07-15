load_data <- function(){readRDS(file = here::here("Data/summarized_major_axis_distribution_data.RData"))}
                      
process_data <- function(data) {
data <- data %>% filter(predicted_species != "none") %>% 
    dplyr::select(-Area, -Major) %>% 
    mutate(day = as.numeric(ymd(date)-ymd("2014-10-06")),
                    fday = as.factor(day),
                    day_range = cut(day, breaks = 6, 
                                     labels = c("0-10d", "11-20d", "21-30d", "31-40d", "41-50d", "50-60d"), right = FALSE),
                    competition = as.factor(ifelse(richness > 1, "poly", "mono")),
                    temp_center = as.vector(scale(temperature, center = T, scale = T))) %>%
  rename(Area = mean_area,
         Major = mean_major,
         Minor = mean_minor)


# calculate individual biovolume 
data <- data %>% mutate(bvol = 4/3 * pi * (Minor/2) ^ 2 * (Major/2))

# add full names for species
data <- data %>% mutate(species_label =
                      case_when(predicted_species == "Colp" ~ "Colpidium striatum",
                                predicted_species == "Dexio" ~ "Dexiostoma campylum",
                                predicted_species == "Loxo" ~ "Loxocephalus sp.",
                                predicted_species == "Para" ~ "Paramecium caudatum",
                                predicted_species == "Spiro" ~ "Spirostomum teres",
                                predicted_species == "Tetra" ~ "Tetrahymena thermophila"))

return(data)                      
}

plot_processed_data <- function(data, rich=1){
  #pdf(here::here(paste0("output/figures/not_for_paper/size_abundance_dynamics_rich_",rich,".pdf")), height=12,width=12)
  p1 <- data %>% group_by(predicted_species, species_label, day, temperature, richness, competition, replicate, combination) %>%  
  dplyr::summarize(mean_size = mean(Major)) %>%
  filter(richness==rich) %>%
  ggplot(data=., aes(x=day, y=mean_size, group=temperature, colour=temperature)) + geom_point() + 
    facet_wrap(combination~ predicted_species, scales="free_y") + 
    geom_smooth(se=F)
  print(p1)
   
  p2 <- data %>% group_by(predicted_species, species_label, day, temperature, richness, competition, replicate, combination) %>%  
    dplyr::summarize(abundance = n()) %>%
    filter(richness==rich) %>%
    ggplot(data=., aes(x=day, y=abundance, group=temperature, colour=temperature)) + geom_point() + geom_smooth(se=F) + 
    facet_wrap(combination~ predicted_species, scales="free_y") + scale_y_log10()
 print(p2)
 #dev.off()
}


create_dataset_analysis_mean_size <- function(data){
 
dd_mean_size_poly <- data %>% 
  group_by(predicted_species, species_label, temperature, richness, competition, replicate, combination, microcosmID) %>%  dplyr::summarize(mean_size = mean(Major, na.rm=T))
  dd_mean_size_poly$sp_rep <- paste0(dd_mean_size_poly$replicate, "_", dd_mean_size_poly$predicted_species)
  
  dd_mean_size_poly$temperature_sq <- dd_mean_size_poly$temperature ^2
  
  dd_mean_size_poly$richness_fac <- factor(dd_mean_size_poly$richness, ordered = F)
  dd_mean_size_poly$temperature_fac <- factor(dd_mean_size_poly$temperature, ordered = F)
  
  dd_mean_size_poly$richness_sq <-dd_mean_size_poly$richness^2
  
  dd_mean_size_poly$rich_sc <- as.numeric(scale(dd_mean_size_poly$richness, center = TRUE, scale = FALSE))
  dd_mean_size_poly$rich_sc_sq <- as.numeric(scale(dd_mean_size_poly$rich_sc^2, center = TRUE, scale = FALSE))
  
  dd_mean_size_poly$temp_sc <- as.numeric(scale(dd_mean_size_poly$temperature, center = TRUE, scale = FALSE))
  dd_mean_size_poly$temp_sc_sq <- as.numeric(scale(dd_mean_size_poly$temp_sc^2, center = TRUE, scale = FALSE))
  
  dd_mean_size_poly <- dd_mean_size_poly %>% mutate(
  incubator = case_when(
    as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
    as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
    as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
    as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
    as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
    as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
    as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
    as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
    as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
    as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
    as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
    as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
  ))
  
  return(dd_mean_size_poly)
}


create_output_folder <- function(){
  
  # check and create directories to store outputs 
  output_dirs <- list("output", "output/tables/", "output/tables/polycultures/",
                      "output/tables/monocultures/", "output/figures/" , "output/figures/main/")

  lapply(output_dirs, function(x) if (!dir.exists(x)) {dir.create(x)})
  
}

twoway_mixed_model_analysis_size <- function(data){

  
  #data <- mean_size_data
  
  nd <- vector("list", 6)
  tw_estimates_lin <- vector("list", 6)
  tw_estimates_poly <- vector("list", 6)
  models <- vector("list", 6)
  aic_select <- vector("list", 6)
  
  species <- c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")
  
  for (i in 1:6){
  
    select_df <- subset(data, species_label == species[i])
    
    mod_length_lin <- lmer(log(mean_size) ~ temp_sc  * rich_sc + (1  | combination) + (1 | incubator),  data=select_df,  REML=F)
    mod_length_lin <- update(mod_length_lin,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    
    # if(species[i] == "Dexiostoma camplyum"){
    #   mod_length_lin <- lmer(log(mean_size) ~ temp_sc  * rich_sc + (1 | combination) + (1|incubator),  data=select_df,  REML=F)
    #   mod_length_lin <- update(mod_length_lin,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    # }
    
    mod_length_poly <- lmer(log(mean_size) ~ (temp_sc_sq + temp_sc) * rich_sc  + (1  | combination) + (1 | incubator),  data=select_df, REML=F)
    mod_length_poly <- update(mod_length_poly,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    
    # if(species[i] == "Dexiostoma camplyum"){
    #   mod_length_poly <- lmer(log(mean_size) ~ (temp_sc_sq + temp_sc) * rich_sc  + (1 | combination) + (1|incubator),  data=select_df, REML=F)
    #   mod_length_poly <- update(mod_length_poly,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    # }
    
    aic_table <- AICcmodavg::aictab(list(mod_length_lin, mod_length_poly), c(
      "y ~ temp", "y ~ temp + temp^2"))
    # Convert to data frame
    aic_select[[i]] <- data.frame(species =species[i], aic_table)
    
    models[[i]] <- list()
    models[[i]] [['Linear mixed model']] <- mod_length_lin
    models[[i]] [['Quadratic mixed model']] <- mod_length_poly
     
    # mods <- list()
    # mods[['Linear mixed model']] <- mod_length_lin
    # mods[['Quadratic mixed model']] <- mod_length_poly
    # msummary(mods, statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
    #         output = here::here("output/tables", paste0(species[i],"_2way_mixed_model_table.html"))) 
     
   
    tw_estimates_lin[[i]] <- broom::tidy(mod_length_lin, conf.int = T)
    tw_estimates_poly[[i]] <- broom::tidy(mod_length_poly, conf.int = T)
    
    nd[[i]] <- ungroup(select_df)  %>% 
      dplyr::select(predicted_species, species_label, combination, competition, temp_sc, temp_sc_sq, rich_sc, rich_sc_sq, richness_fac, mean_size, incubator) %>% 
      group_by(predicted_species,species_label, combination, competition, temp_sc, temp_sc_sq,rich_sc, rich_sc_sq, richness_fac, incubator, .add=F) %>%
      dplyr::summarize(mean_size = mean(mean_size))
    
    nd[[i]]$pred_size <- predict(mod_length_poly, newdata = nd[[i]], re.form = ~(1 | combination))
    
  }

  cm <- c('(Intercept)' = 'Constant',
          'rich_sc' = 'Scaled richness',
          'rich_sc_sq' = 'Scaled richness squared',       
          'temp_sc' = 'Scaled temperature',
          'temp_sc_sq' = 'Scaled temperature squared',
          'temp_sc:rich_sc' = 'Temperature : richness',
          'temp_sc_sq:rich_sc' = 'Temperature squared: richness',
          'cor__(Intercept).temp_sc' = "Intercept-slope correlation",
          'sd__(Intercept)' = 'SD Intercept',
          'sd__Observation' = "SD observation",
          'sd__temp_sc' = "SD temperature slope"
  )
  # give names to models
  names(models) <- species
  
  #  msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #           output = here::here("output/tables/polycultures/2way_mixed_model_table_size_linear.docx"))
  #  
  # msummary(lapply(models, "[[", 2), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #          output = here::here("output/tables/polycultures/2way_mixed_model_table_size_quadratic.docx"))
  
  model_obj <- c(lapply(models, "[[", 2)[c(3,6)], lapply(models, "[[", 1)[-c(3,6)])
  model_obj <- model_obj[c(3,4,1,5,6,2)]
  msummary(model_obj, statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
           output = here::here("output/tables/polycultures/2way_mixed_model_table_size.docx"))
  
  aic_select_df <- bind_rows(aic_select)
  
  highlight_rows <- aic_select_df %>%
    group_by(species) %>%
    # Filter to only those with Delta_AICc < 2
    filter(Delta_AICc < 2) %>%
    # For each set, find the model with minimum K
    slice_min(K) %>%
    # Create a unique identifier for each row
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create the highlight vector
  rows_to_highlight <- highlight_rows$row_id
  
  # Add a temporary identifier column to the main data
  all_tables <- aic_select_df %>%
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create a gt table with grouped headers
  final_table <- all_tables %>%
    gt() %>%
    tab_header(
      title = "Model Comparison across Multiple Datasets",
      subtitle = "AICc Model Selection Results"
    ) %>%
    tab_options(
      row_group.background.color = "#f7f7f7",
      heading.background.color = "#e8e8e8",
      column_labels.background.color = "#e0e0e0"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) %>%
    # Add bold formatting to rows that meet our criteria
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = row_id %in% rows_to_highlight
      )
    ) %>%
    fmt_number(
      columns = c("AICc", "Delta_AICc", "AICcWt", "Cum.Wt", "LL"),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c("ModelLik"),
      decimals = 3
    ) %>%
    cols_align(
      align = "center",
      columns = -c("Modnames")
    ) %>%
    # Remove the temporary row_id column
    cols_hide(columns = "row_id")
  # To export to Word
  # Option 1: Using gt and webshot2 approach
  gt::gtsave(final_table, "output/tables/polycultures/aic_tables_size.docx")
  
  
  
  nd_df <- bind_rows(nd)
  
  # check summary
  #lapply(1:length(models), function(x) summary(models[[x]][1]$`Linear mixed model`))
  #lapply(1:length(models), function(x) summary(models[[x]][2]$`Quadratic mixed model`))

  
  # bootstrapping to get confidence intervals of model parameter (including random effects)
  #confint(models[[1]][1]$`Linear mixed model`, method="boot", nsim=500, oldNames= FALSE)
  bootstrap_confidence_list <- lapply(1:6, function(x) confint(models[[x]][2]$`Quadratic mixed model`, method="boot", nsim=500, oldNames= FALSE))
  bootstrap_confidence_df <- lapply(bootstrap_confidence_list, as.data.frame)
  bootstrap_confidence_df <- bind_rows(bootstrap_confidence_df)

  terms <- c("sd__(Intercept combination)", "sd__(Intercept incubator)", "sd__Observation", "(Intercept)", "temp_sc_sq", "temp_sc", "rich_sc", "temp_sc_sq:rich_sc", "temp_sc:rich_sc")
  bootstrap_confidence_df$term <- rep(terms, times=6)
  bootstrap_confidence_df$species_label <- rep(species, each=length(terms))
  
  tw_estimates_df <- bind_rows(tw_estimates_poly) %>% 
    mutate(species_label = rep(species, times = length(terms) ))
  
  
  
  
  tw_estimates_df <- merge(tw_estimates_df, bootstrap_confidence_df, by=c("term", "species_label"))
  tw_estimates_df <- tw_estimates_df %>% mutate(
    term = case_when(
      term == "(Intercept)" ~ "constant",
      term == 'temp_sc' ~ 'scaled temperature',
      term == 'rich_sc' ~ 'scaled richness',
      term == 'temp_sc_sq' ~ 'scaled temperature squared',
      term == 'temp_sc:rich_sc' ~ 'scaled temperature x richness',
      term == 'temp_sc_sq:rich_sc' ~ 'scaled temperature squared x richness',
      term == 'sd__(Intercept)' ~ 'SD intercept (combination)',
      term == 'sd__temp_sc' ~ 'SD slope',
      term == 'cor__(Intercept combination).temp_sc' ~ 'Slope intercept correlation',
      term == 'sd__(Intercept incubator)' ~ 'SD intercept (incubator)',
      term == 'sd__Observation' ~ 'SD observation.level',
      TRUE ~ term,
    )
  )
  return(list(tw_estimates_df, nd_df, models))
}


linear_model_analysis_size <- function(data){
  
  #data <- mean_size_data
  
  nd <- vector("list", 6)
  tw_estimates_lin <- vector("list", 6)
  tw_estimates_poly <- vector("list", 6)
  models <- vector("list", 6)
  aic_select <- vector("list", 6)
  
  species <- c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")

  #pdf("output/figures/not_for_paper/model_diagnostics_size.pdf")
  
  for (i in 1:6){
    
    select_df <- subset(data, species_label == species[i] & richness == 1) 
    # uncomment to aggregate data points per incubator
    #%>% group_by(incubator, temp_sc, temp_sc_sq) %>% summarize(mean_size = mean(mean_size, na.rm=T))
    
    # incubator RE
    #mod_length_lin <- lmer(log(mean_size) ~ temp_sc + (1|incubator), REML=F, data=select_df)
    #mod_length_poly <- lmer(log(mean_size) ~ (temp_sc_sq + temp_sc) + (1|incubator), REML=F, data=select_df)
    
    mod_length_lin <- lm(log(mean_size) ~ temp_sc, data=select_df)
    mod_length_poly <- lm(log(mean_size) ~ (temp_sc + temp_sc_sq), data=select_df)
    
    aic_table <- AICcmodavg::aictab(list(mod_length_lin, mod_length_poly), c(
      "y ~ temp", "y ~ temp + temp^2"))
    # Convert to data frame
    aic_select[[i]] <- data.frame(species =species[i], aic_table)

    
    models[[i]] <- list()
    models[[i]] [['Linear mixed model']] <- mod_length_lin
    models[[i]] [['Quadratic mixed model']] <- mod_length_poly
    
  }


  #dev.off()
  # cm <- c('(Intercept)' = 'Constant',
  #         'temp_sc' = 'Scaled temperature',
  #         'temp_sc_sq' = 'Scaled temperature squared',
  #         'sd__observation' = 'SD intercept')
 
   # give names to models
  names(models) <- species
  
   # msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
   #          output = here::here("output/tables/monocultures/linear_model_table_size_linear.docx"))
   # 
   # msummary(lapply(models, "[[", 2), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
   #         output = here::here("output/tables/monocultures/linear_model_table_size_quadratic.docx"))
  
    msummary(c(lapply(models, "[[", 2)[1:2], lapply(models, "[[", 1)[-c(1:2)]), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
             output = here::here("output/tables/monocultures/linear_model_table_size.docx"))
    
    aic_select_df <- bind_rows(aic_select)
    
    highlight_rows <- aic_select_df %>%
      group_by(species) %>%
      # Filter to only those with Delta_AICc < 2
      filter(Delta_AICc < 2) %>%
      # For each set, find the model with minimum K
      slice_min(K) %>%
      # Create a unique identifier for each row
      mutate(row_id = paste(species, Modnames, sep = "_"))
    
    # Create the highlight vector
    rows_to_highlight <- highlight_rows$row_id
    
    # Add a temporary identifier column to the main data
    all_tables <- aic_select_df %>%
      mutate(row_id = paste(species, Modnames, sep = "_"))
    
    # Create a gt table with grouped headers
    final_table <- all_tables %>%
      gt() %>%
      tab_header(
        title = "Model Comparison across Multiple Datasets",
        subtitle = "AICc Model Selection Results"
      ) %>%
      tab_options(
        row_group.background.color = "#f7f7f7",
        heading.background.color = "#e8e8e8",
        column_labels.background.color = "#e0e0e0"
      ) %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_row_groups()
      ) %>%
      # Add bold formatting to rows that meet our criteria
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(
          rows = row_id %in% rows_to_highlight
        )
      ) %>%
      fmt_number(
        columns = c("AICc", "Delta_AICc", "AICcWt", "Cum.Wt", "LL"),
        decimals = 2
      ) %>%
      fmt_number(
        columns = c("ModelLik"),
        decimals = 3
      ) %>%
      cols_align(
        align = "center",
        columns = -c("Modnames")
      ) %>%
      # Remove the temporary row_id column
      cols_hide(columns = "row_id")
    # To export to Word
    # Option 1: Using gt and webshot2 approach
    gt::gtsave(final_table, "output/tables/monocultures/aic_tables_size.docx")
    
   
  # check summary
  #lapply(1:length(models), function(x) summary(models[[x]][1]$`Linear mixed model`))
  #lapply(1:length(models), function(x) summary(models[[x]][2]$`Quadratic mixed model`))
  
}


linear_model_analysis_supply <- function(data, supply_proxy = "max_bm"){
  #browser()
  
  
  data <- data %>% mutate(richness_sq = richness^2,
                          rich_sc = as.numeric(scale(richness, center = TRUE, scale = FALSE)),
                          rich_sc_sq = as.numeric(scale(rich_sc^2, center = TRUE, scale = FALSE)),
                          temp_sq = richness^2,
                          temp_sc = as.numeric(scale(temperature, center = TRUE, scale = FALSE)),
                          temp_sc_sq = as.numeric(scale(temp_sc^2, center = TRUE, scale = FALSE)),
                          incubator = case_when(
                            as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
                            as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
                            as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
                            as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
                            as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
                            as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
                            as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
                            as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
                            as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
                            as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
                            as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
                            as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
                          ))
  nd <- vector("list", 6)
  tw_estimates_lin <- vector("list", 6)
  tw_estimates_poly <- vector("list", 6)
  models <- vector("list", 6)
  aic_select <- vector("list", 6)
  
  species <- c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")
  
  #pdf("output/figures/not_for_paper/model_diagnostics_supply.pdf")
  for (i in 1:6){
    
    select_df <- subset(data, species_label == species[i] & richness == 1)
    
    # mod_length_lin <- lmer(log(get(supply_proxy)) ~ temp_sc + (1|incubator), REML=F,  data=select_df)
    # mod_length_poly <- lmer(log(get(supply_proxy)) ~ (temp_sc_sq + temp_sc) + (1|incubator),  REML=F, data=select_df)
    # print(AICcmodavg::aictab(list(mod_length_lin, mod_length_poly)))
    
    mod_length_lin <- lm(log(get(supply_proxy)) ~ temp_sc,  data=select_df)
    mod_length_poly <- lm(log(get(supply_proxy)) ~ (temp_sc + temp_sc_sq), data=select_df)
  
    aic_table <- AICcmodavg::aictab(list(mod_length_lin, mod_length_poly), c(
      "y ~ temp", "y ~ temp + temp^2"))
    aic_select[[i]] <- data.frame(species =species[i], aic_table)
    
    models[[i]] <- list()
    models[[i]] [['Linear mixed model']] <- mod_length_lin
    models[[i]] [['Quadratic mixed model']] <- mod_length_poly
    
  }
    #dev.off()
  
  cm <- c('(Intercept)' = 'Constant',
          'temp_sc' = 'Scaled temperature',
          'temp_sc_sq' = 'Scaled temperature squared')
  # give names to models
  names(models) <- species
  
  #  msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #           output = here::here("output/tables/monocultures/linear_model_table_supply_linear.docx"))
  #  
  # msummary(lapply(models, "[[", 2), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #          output = here::here("output/tables/monocultures/linear_model_table_supply_quadratic.docx"))
  #browser()
  msummary(c(lapply(models, "[[", 2)[1], lapply(models, "[[", 1)[-1]), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
           output = here::here("output/tables/monocultures/linear_model_table_supply.docx"))
  
  aic_select_df <- bind_rows(aic_select)
  
  highlight_rows <- aic_select_df %>%
    group_by(species) %>%
    # Filter to only those with Delta_AICc < 2
    filter(Delta_AICc < 2) %>%
    # For each set, find the model with minimum K
    slice_min(K) %>%
    # Create a unique identifier for each row
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create the highlight vector
  rows_to_highlight <- highlight_rows$row_id
  
  # Add a temporary identifier column to the main data
  all_tables <- aic_select_df %>%
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create a gt table with grouped headers
  final_table <- all_tables %>%
    gt() %>%
    tab_header(
      title = "Model Comparison across Multiple Datasets",
      subtitle = "AICc Model Selection Results"
    ) %>%
    tab_options(
      row_group.background.color = "#f7f7f7",
      heading.background.color = "#e8e8e8",
      column_labels.background.color = "#e0e0e0"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) %>%
    # Add bold formatting to rows that meet our criteria
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = row_id %in% rows_to_highlight
      )
    ) %>%
    fmt_number(
      columns = c("AICc", "Delta_AICc", "AICcWt", "Cum.Wt", "LL"),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c("ModelLik"),
      decimals = 3
    ) %>%
    cols_align(
      align = "center",
      columns = -c("Modnames")
    ) %>%
    # Remove the temporary row_id column
    cols_hide(columns = "row_id")
  # To export to Word
  # Option 1: Using gt and webshot2 approach  
  gt::gtsave(final_table, "output/tables/monocultures/aic_tables_supply.docx")
  
  
  return(models)
  
}


linear_model_analysis_demand <- function(data){
  
  #browser()
  #data <- logistic_growth_parameters_as_demand_proxy
  
  data <- data %>% filter(term == "r") %>% mutate(incubator = case_when(
                            as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
                            as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
                            as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
                            as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
                            as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
                            as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
                            as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
                            as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
                            as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
                            as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
                            as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
                            as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
                          ))
  
  nd <- vector("list", 6)
  tw_estimates_lin <- vector("list", 6)
  tw_estimates_poly <- vector("list", 6)
  models <- vector("list", 6)
  aic_select <- vector("list", 6)
  
  species <- c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")
  
  for (i in 1:6){
    
    select_df <- subset(data, species_label == species[i]) 
    #%>% group_by(temp_sc, temp_sc_sq, incubator) %>% mutate(estimate = mean(estimate, na.rm=T))

    
    # mod_length_lin <- lmer(log(estimate) ~ temp_sc + (1|incubator),  data=select_df, REML=F)
    # mod_length_poly <- lmer(log(estimate) ~ temp_sc + temp_sc_sq + (1|incubator),  data=select_df, REML=F)
    # print(AICcmodavg::aictab(list(mod_length_lin, mod_length_poly)))
    
    mod_length_lin <- lm(log(estimate) ~ temp_sc,  data=select_df)
    mod_length_poly <- lm(log(estimate) ~ temp_sc + temp_sc_sq,  data=select_df)
  
      aic_table <- AICcmodavg::aictab(list(mod_length_lin, mod_length_poly), c(
      "y ~ temp", "y ~ temp + temp^2"))
    # Convert to data frame
    aic_select[[i]] <- data.frame(species =species[i], aic_table)
    
    models[[i]] <- list()
    models[[i]] [['Linear mixed model']] <- mod_length_lin
    models[[i]] [['Quadratic mixed model']] <- mod_length_poly
    
  }
  
  
  cm <- c('(Intercept)' = 'Constant',
          'temp_sc' = 'Scaled temperature',
          'temp_sc_sq' = 'Scaled temperature squared'  )
  
  # give names to models
  names(models) <- species
  
  # msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #           output = here::here("output/tables/monocultures/linear_model_table_demand_linear.docx"))
  #  
  # msummary(lapply(models, "[[", 2), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #          output = here::here("output/tables/monocultures/linear_model_table_demand_quadratic.docx"))

  msummary(c(lapply(models, "[[", 2)[1], lapply(models, "[[", 1)[-1]), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
           output = here::here("output/tables/monocultures/linear_model_table_demand.docx"))
  
  aic_select_df <- bind_rows(aic_select)
  
  highlight_rows <- aic_select_df %>%
    group_by(species) %>%
    # Filter to only those with Delta_AICc < 2
    filter(Delta_AICc < 2) %>%
    # For each set, find the model with minimum K
    slice_min(K) %>%
    # Create a unique identifier for each row
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create the highlight vector
  rows_to_highlight <- highlight_rows$row_id
  
  # Add a temporary identifier column to the main data
  all_tables <- aic_select_df %>%
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create a gt table with grouped headers
  final_table <- all_tables %>%
    gt() %>%
    tab_header(
      title = "Model Comparison across Multiple Datasets",
      subtitle = "AICc Model Selection Results"
    ) %>%
    tab_options(
      row_group.background.color = "#f7f7f7",
      heading.background.color = "#e8e8e8",
      column_labels.background.color = "#e0e0e0"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) %>%
    # Add bold formatting to rows that meet our criteria
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = row_id %in% rows_to_highlight
      )
    ) %>%
    fmt_number(
      columns = c("AICc", "Delta_AICc", "AICcWt", "Cum.Wt", "LL"),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c("ModelLik"),
      decimals = 3
    ) %>%
    cols_align(
      align = "center",
      columns = -c("Modnames")
    ) %>%
    # Remove the temporary row_id column
    cols_hide(columns = "row_id")
  # To export to Word
  # Option 1: Using gt and webshot2 approach
  gt::gtsave(final_table, "output/tables/monocultures/aic_tables_demand.docx")
  
  
}


twoway_mixed_model_analysis_SD_ratio <- function(data){
  
  #data <- mean_size_data

  data <- data %>% group_by(predicted_species) %>% 
    mutate(log_supply_over_log_demand_c = log_supply_over_log_demand-mean(log_supply_over_log_demand),
           rich_sc = richness - mean(richness, na.rm=T), rich_sc_sq = rich_sc^2,  incubator = case_when(
             as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
             as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
             as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
             as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
             as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
             as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
             as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
             as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
             as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
             as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
             as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
             as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
           ))
  data$rich_fac <- as.factor(data$richness)
  
  nd <- vector("list", 6)
  tw_estimates_lin <- vector("list", 6)
  tw_estimates_poly <- vector("list", 6)
  models <- vector("list", 6)
  aic_select <- vector("list", 6)
  
  species <- c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")
  
  for (i in 1:6){
    
    select_df <- subset(data, species_label == species[i])
    
    mod_length_lin <- lmer(mean_log_major ~ log_supply_over_log_demand_c * rich_sc + (1 | combination) + (1|incubator), 
                 data=select_df, REML=F)
    mod_length_lin <- update(mod_length_lin,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    
    mod_length_poly <- lmer(mean_log_major ~ log_supply_over_log_demand_c * rich_sc + rich_sc_sq + (1 | combination) + (1|incubator),  data=select_df, REML=F)
 

    models[[i]] <- list()
    models[[i]] [['Linear mixed model']] <- mod_length_lin
    models[[i]] [['Quadratic mixed model']] <- mod_length_poly
    
    # mods <- list()
    # mods[['Linear mixed model']] <- mod_length_lin
    # mods[['Quadratic mixed model']] <- mod_length_poly
    # msummary(mods, statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
    #         output = here::here("output/tables", paste0(species[i],"_2way_mixed_model_table.html"))) 
    
    
    tw_estimates_lin[[i]] <- broom::tidy(mod_length_lin, conf.int = T)
    tw_estimates_poly[[i]] <- broom::tidy(mod_length_poly, conf.int = T)
    
    #browser()
    
    aic_table <- AICcmodavg::aictab(list(mod_length_lin, mod_length_poly), c(
      "y ~ rich", "y ~ rich + rich^2"))
    # Convert to data frame
    aic_select[[i]] <- data.frame(species =species[i], aic_table)
    
    nd[[i]] <- ungroup(select_df)  %>% 
      dplyr::select(predicted_species, species_label, combination, rich_sc, rich_sc_sq, mean_log_major, log_supply_over_log_demand_c) %>% 
      group_by(predicted_species,species_label, combination, rich_sc, rich_sc_sq, log_supply_over_log_demand_c, .add=F) %>%
      dplyr::summarize(mean_log_major = mean(mean_log_major))
    
    nd[[i]]$pred_mean_log_major <- predict(mod_length_poly, newdata = nd[[i]], re.form=NA)
    
  }
  
  cm <- c('(Intercept)' = 'Constant',
          'log_supply_over_log_demand_c' = 'log supply over log demand',
          'rich_sc' = 'Scaled richness',
          'log_supply_over_log_demand_c:rich_sc' = 'log supply over log demand x scaled richness squared',
          'rich_sc_sq' = 'Scaled richness squared',       
          'sd__(Intercept)' = 'SD Intercept',
          'sd__Observation' = "SD observation",
          'sd__rich_sc' = "SD richness slope"
  )
  # give names to models
  names(models) <- species
  
  #browser()
  
   msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
            output = here::here("output/tables/polycultures/2way_mixed_model_table_SD.docx"))
   
  # msummary(lapply(models, "[[", 2), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
  #          output = here::here("output/tables/polycultures/2way_mixed_model_table_SD_quadratic.docx"))
  # 
  
  aic_select_df <- bind_rows(aic_select)
  
  highlight_rows <- aic_select_df %>%
    group_by(species) %>%
    # Filter to only those with Delta_AICc < 2
    filter(Delta_AICc < 2) %>%
    # For each set, find the model with minimum K
    slice_min(K) %>%
    # Create a unique identifier for each row
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create the highlight vector
  rows_to_highlight <- highlight_rows$row_id
  
  # Add a temporary identifier column to the main data
  all_tables <- aic_select_df %>%
    mutate(row_id = paste(species, Modnames, sep = "_"))
  
  # Create a gt table with grouped headers
  final_table <- all_tables %>%
    gt() %>%
    tab_header(
      title = "Model Comparison across Multiple Datasets",
      subtitle = "AICc Model Selection Results"
    ) %>%
    tab_options(
      row_group.background.color = "#f7f7f7",
      heading.background.color = "#e8e8e8",
      column_labels.background.color = "#e0e0e0"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) %>%
    # Add bold formatting to rows that meet our criteria
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = row_id %in% rows_to_highlight
      )
    ) %>%
    fmt_number(
      columns = c("AICc", "Delta_AICc", "AICcWt", "Cum.Wt", "LL"),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c("ModelLik"),
      decimals = 3
    ) %>%
    cols_align(
      align = "center",
      columns = -c("Modnames")
    ) %>%
    # Remove the temporary row_id column
    cols_hide(columns = "row_id")
  # To export to Word
  # Option 1: Using gt and webshot2 approach
  gt::gtsave(final_table, "output/tables/polycultures/aic_tables_SD_ratio.docx")
  
  
  nd_df <- bind_rows(nd)
  
  # check summary
  #lapply(1:length(models), function(x) summary(models[[x]][1]$`Linear mixed model`))
  #lapply(1:length(models), function(x) summary(models[[x]][2]$`Quadratic mixed model`))
  #browser()
  tw_estimates_df <- bind_rows(tw_estimates_poly) %>% 
    mutate(species_label = rep(species, times = c(8,8,8,8,8,8)))
  
  
  # bootstrapping to get confidence intervals of model parameter (including random effects)
  #confint(models[[1]][1]$`Linear mixed model`, method="boot", nsim=500, oldNames= FALSE)
  bootstrap_confidence_list <- lapply(1:6, function(x) confint(models[[x]][2]$`Quadratic mixed model`, method="boot", nsim=500, oldNames= FALSE))
  bootstrap_confidence_df <- lapply(bootstrap_confidence_list, as.data.frame)
  bootstrap_confidence_df <- bind_rows(bootstrap_confidence_df)
  
  terms <- c("sd__(Intercept combination)", "sd__(Intercept incubator)", "sd__Observation", "(Intercept)", "log_supply_over_log_demand_c","rich_sc", "rich_sc_sq", "log_supply_over_log_demand_c:rich_sc")
  
  bootstrap_confidence_df$term <- rep(terms, times=6)
  bootstrap_confidence_df$species_label <- rep(species, each=8)
  
  
  tw_estimates_df <- merge(tw_estimates_df, bootstrap_confidence_df, by=c("term", "species_label"))
  
  tw_estimates_df <- tw_estimates_df %>% mutate(
    term = case_when(
      term == "(Intercept)" ~ "constant",
      term == 'rich_sc' ~ 'scaled richness',
      term == 'rich_sc_sq' ~ 'scaled richness squared',
      term == 'log_supply_over_log_demand_c' ~ 'log supply over log demand scaled',
      term == 'log_supply_over_log_demand_c:rich_sc' ~ 'log supply over log demand x scaled richness',
      term == 'sd__(Intercept)' ~ 'SD intercept',
      term == 'sd__Observation' ~ 'SD observation.level',
      TRUE ~ term,
    )
  )
  return(list(tw_estimates_df, nd_df, models))
}



fit_logistic_growth <- function(data, response_var){

  pop_dyn <- data %>% group_by(predicted_species, temperature, day, replicate, microcosmID) %>%
    filter(richness == 1) %>% dplyr::summarize(count = n())
  
  
  microcosmID <- data %>% group_by(predicted_species, temperature, replicate, microcosmID) %>%
    filter(richness == 1) %>% dplyr::summarize(microcosmID = unique(microcosmID))

  initial_conditions <- expand.grid(predicted_species=unique(pop_dyn$predicted_species), 
                                    replicate=unique(pop_dyn$replicate),
                                    temperature=unique(pop_dyn$temperature))

  # merge initial conditions with microcosmID to preserve microcosmID info
  initial_conditions <- merge(microcosmID, initial_conditions)
  initial_conditions$day <- 0
  initial_conditions$response_var <- 0
  
  if(response_var == "density"){
    # set initial density to 3 and scale up to individuals per ml
    initial_conditions$response_var <- 3/100/23
    
    mono_dyn <- data %>% group_by(predicted_species, temperature, day, replicate, microcosmID) %>%
      filter(richness == 1) %>% dplyr::summarize(response_var = n())
    
  } else {
    bvol_df <- data %>% filter(richness == 1) %>% group_by(predicted_species, temperature) %>% dplyr::summarize(mean_bvol = mean(bvol, na.rm = T))
    initial_conditions <- merge(initial_conditions, bvol_df)
    initial_conditions$response_var <- initial_conditions$mean_bvol * 3/100/23
    initial_conditions$mean_bvol <- NULL
    
    mono_dyn <- data %>% group_by(predicted_species, temperature, day, replicate) %>%
      filter(richness == 1) %>% dplyr::summarize(response_var = sum(bvol, na.rm = T))
    
  }
  
  mono_dyn <- rbind(as.data.frame(mono_dyn), initial_conditions)
  
  mono_dyn_ts_cut <- mono_dyn %>% mutate(keep = case_when(predicted_species == "Colp" & day < 20 ~ 1,
                                                          predicted_species == "Dexio" & day < 20 ~ 1,
                                                          predicted_species == "Loxo" & day < 25 ~ 1,
                                                          predicted_species == "Para" & day < 20 ~ 1,
                                                          predicted_species == "Spiro" & day < 40 ~ 1,
                                                          predicted_species == "Tetra" & day < 10 ~ 1)) %>%
    filter(keep == 1) %>%
    arrange(predicted_species, temperature, day, replicate, microcosmID)
  
  
  fits <- mono_dyn_ts_cut %>%
    group_by(predicted_species, temperature, replicate, microcosmID) %>%
    nest() %>%
    mutate(fit = purrr::map(data, ~ nls_multstart(response_var ~ response_var[1] * K / (response_var[1] + ((K - response_var[1]) * exp(-r * day))),
                                                  data = .x,
                                                  iter = 500,
                                                  start_lower = c(K = 1, r = 0.1),
                                                  start_upper = c(K = 1000, r = 10),
                                                  supp_errors = 'Y',
                                                  na.action = na.omit,
                                                  lower = c(K= 0.1, r = 0.0001))))
  
  
  # quick out estimates with improper input parameters values
  #fits2 <- fits[-c(9, 20, 45, 59, 62, 74, 79, 89,97),]
  fits2 <- fits[-c(20, 57),]
  
  # this is needed because of code breaking changes to nest / unnest
  nest <- nest_legacy 
  unnest <- unnest_legacy
  
  # get summary info
  info <- fits2 %>% unnest(fit %>% map(glance))
  
  # get params
  params <- fits2 %>% unnest(fit %>% map(tidy))
  
  # get confidence intervals
  CI <- fits2 %>% 
    unnest(fit %>% map(~ confint2(.x) %>%
                         data.frame() %>%
                         rename(., conf.low = X2.5.., conf.high = X97.5..))) %>%
    group_by(., temperature, predicted_species, replicate, microcosmID) %>%
    mutate(., term = c('r', 'K')) %>%
    ungroup()

  # merge parameters and CI estimates
  params <- merge(params, CI, by = intersect(names(params), names(CI)))

  
  # get predictions
  preds <- fits2 %>%
    unnest(fit %>% map(augment))
  
  new_preds <- mono_dyn_ts_cut %>%
    do(data.frame(day = seq(1, 40, length.out = 150), stringsAsFactors = FALSE)) 
  
  # create new predictions
  preds2 <- fits2 %>%
    unnest(fit %>% map(augment, newdata = new_preds)) %>%
    rename(response_var = .fitted) %>%
    ungroup() %>% 
    mutate(keep = case_when(predicted_species == "Colp" & day < 20 ~ 1,
                            predicted_species == "Dexio" & day < 20 ~ 1,
                            predicted_species == "Loxo" & day < 25 ~ 1,
                            predicted_species == "Para" & day < 20 ~ 1,
                            predicted_species == "Spiro" & day < 40 ~ 1,
                            predicted_species == "Tetra" & day < 10 ~ 1)) %>%
    filter(keep == 1) %>% dplyr::select(-keep)
  
  # add full names for species
  params <- params %>% mutate(species_label =
                                          case_when(predicted_species == "Colp" ~ "Colpidium striatum",
                                                    predicted_species == "Dexio" ~ "Dexiostoma campylum",
                                                    predicted_species == "Loxo" ~ "Loxocephalus sp.",
                                                    predicted_species == "Para" ~ "Paramecium caudatum",
                                                    predicted_species == "Spiro" ~ "Spirostomum teres",
                                                    predicted_species == "Tetra" ~ "Tetrahymena thermophila")
  )
  
  # add scaled and squared terms
  params$temp_sc <- as.numeric(scale(params$temperature, center=T, scale=F))
  params$temp_sc_sq <- as.numeric(scale(params$temp_sc^2, center=T, scale=F))
  
  return(params)
}

align_mono_and_polyculture_time_series <- function(data, params, response_var){

  # first identify the day where the monocultures cross the  threshold
  pop_dyn <- data %>% group_by(predicted_species, temperature, day, replicate) %>%
    filter(richness == 1) %>% dplyr::summarize(count = n())
  
  mono_K <- params %>% filter(term == "K") %>% dplyr::select(predicted_species, temperature, estimate)
  mono_K <- merge(mono_K, data.frame(replicate=1:3))
  mono_K <- mono_K %>% mutate(dens_20 = estimate/5)
  
  pop_dyn_cut <- merge(pop_dyn, mono_K) 
  pop_dyn_cut <- pop_dyn_cut %>% group_by(predicted_species, temperature, replicate) %>% mutate(keep = ifelse(count > dens_20, 1, 0))
  
  keep_days <- pop_dyn_cut %>% group_by(predicted_species, temperature, replicate) %>% filter(keep==1) %>% top_n(day, n=-1) %>% rename(keep_day = day) %>% dplyr::select(predicted_species, temperature, replicate, keep_day)

initial_conditions <- expand.grid(predicted_species=unique(pop_dyn$predicted_species), 
                                    replicate=unique(pop_dyn$replicate),
                                    temperature=unique(pop_dyn$temperature))
  initial_conditions$day <- 0
  
  if(response_var == "density"){
    initial_conditions$response_var <- 3/100/23
    
    mono_dyn <- data %>% group_by(predicted_species, temperature, day, replicate) %>%
      filter(richness == 1) %>% dplyr::summarize(response_var = n())
    
  } else {
    bvol_df <- data %>% filter(richness == 1) %>% group_by(predicted_species, temperature) %>% dplyr::summarize(mean_bvol = mean(bvol, na.rm = T))
    initial_conditions <- merge(initial_conditions, bvol_df)
    initial_conditions$response_var <- initial_conditions$mean_bvol * 3/100/23
    initial_conditions$mean_bvol <- NULL
    
    mono_dyn <- data %>% group_by(predicted_species, temperature, day, replicate) %>%
      filter(richness == 1) %>% dplyr::summarize(response_var = sum(bvol, na.rm = T))
    
  }
  
  mono_dyn <- rbind(as.data.frame(mono_dyn), initial_conditions)
  
  mono_dyn_cut <- merge(mono_dyn, keep_days)
  mono_dyn_cut <- mono_dyn_cut %>% filter(day>keep_day)
  
  # reset first day of each cut time series
  mono_dyn_cut <- mono_dyn_cut %>% mutate(day = day-keep_day)
  
  return(mono_dyn_cut)
}


# define function to calculate area under the curve
area_under_curve <- function(response_var, time) {
  id <- order(time)
  sum(diff(time[id])*rollmean(response_var[id],2))
}

calculate_supply_proxy_mono <- function(aligned_time_series, processed_data){

# work with aligned time series
dd_auc <- aligned_time_series %>% filter(day <= 40) %>% arrange(predicted_species, temperature, replicate, day) %>% group_by(predicted_species, temperature, replicate) %>% nest()

# calculate AUC
dd_auc$AUC <- unlist(lapply(1:nrow(dd_auc), function(x) area_under_curve(response_var = dd_auc$data[[x]]$response_var, time=dd_auc$data[[x]]$day)))
  
# calculate max biomass
dd_auc$max_bm <- unlist(lapply(1:nrow(dd_auc), function(x) max(dd_auc$data[[x]]$response_var, na.rm=T)))
  
 dd_species_mean_size <- processed_data %>%  filter(richness == 1) %>% 
     group_by(temperature, combination, richness, competition, predicted_species, replicate) %>% 
     dplyr::summarize(mean_log_major = mean(log(Major), na.rm=T))
  
dd_auc <- merge(dd_auc, dd_species_mean_size)

return(dd_auc)
}


create_mono_SD_data <- function(demand_data, supply_data, supply_proxy = "max_bm"){
  
  demand_data <- demand_data %>% filter(term == "r")
  
  dd_auc2 <- merge(supply_data, demand_data, by = c("predicted_species", "temperature", "replicate"))
  
  # uncomment to base results on AUC
  #dd_auc2 <- dd_auc2 %>% filter(AUC > 0) %>% mutate(log_AUC_over_log_r = log(AUC)-log(estimate))
  dd_auc2 <- dd_auc2 %>% mutate(log_supply_over_log_demand = log(get(supply_proxy))-log(estimate))
  
  # add full names for species
  dd_auc2 <- dd_auc2 %>% mutate(species_label =
                                  case_when(predicted_species == "Colp" ~ "Colpidium striatum",
                                            predicted_species == "Dexio" ~ "Dexiostoma campylum",
                                            predicted_species == "Loxo" ~ "Loxocephalus sp.",
                                            predicted_species == "Para" ~ "Paramecium caudatum",
                                            predicted_species == "Spiro" ~ "Spirostomum teres",
                                            predicted_species == "Tetra" ~ "Tetrahymena thermophila"),
                                incubator = case_when(
                                  as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
                                  as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
                                  as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
                                  as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
                                  as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
                                  as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
                                  as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
                                  as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
                                  as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
                                  as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
                                  as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
                                  as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
                                ))
  
  return(dd_auc2)
}



linear_model_analysis_SD <- function(data){ 

mods <- vector("list", 6)
models <- vector("list", 6)
dds <- vector("list", 6)

counter <- 1

for (i in c("Colp", "Dexio", "Loxo", "Para", "Spiro", "Tetra")){

#pdf(paste0("output/figures/not_for_paper/model_diagnostics_SD_",i ,".pdf"))
selected_species <- i

dd <- subset(data, predicted_species == selected_species) %>% mutate(log_supply_over_log_demand_c = log_supply_over_log_demand - mean(log_supply_over_log_demand, na.rm=T))

# mods[[counter]] <- lmer(mean_log_major ~ log_supply_over_log_demand_c + (1|incubator), data=dd)
mods[[counter]] <- lm(mean_log_major ~ log_supply_over_log_demand_c, data=dd)
#print(autoplot(mods[[counter]], which = 1:6, ncol = 3, label.size = 3))

cm <- c('(Intercept)' = 'Constant',
        'log_supply_over_log_demand_c' = 'log(S/D)'
)
#dev.off()

models[[counter]] <- list()
models[[counter]][['Linear model']] <- mods[[counter]]

data2 <- dd
data2$mean_log_major <- predict(mods[[counter]], data2)
dds[[counter]] <- data2

counter <- counter+1


}



#print(lapply(mods, glance))

names(models) <-  c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")

msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC', 
         output = here::here("output/tables/monocultures/linear_model_table_SD.docx"))

return(as.data.frame(bind_rows(dds)))
}




build_composite_figure_SD_model_mono <- function(data){
  
  p0 <- ggplot() + 
    geom_point(data=data, aes(x=temperature, y=mean_log_major, colour=species_label, group=species_label, fill=species_label), size = 1.5)  + 
    facet_grid(species_label~., scale = "free", labeller = label_wrap_gen(width=10)) + 
    ylab("Mean log cell length (in micrometers)") +
    stat_smooth(data=filter(data, species_label %in% c("Colpidium striatum", "Dexiostoma campylum")), aes(x=temperature, y=mean_log_major, colour=species_label, group=species_label, fill=species_label), method = "lm", formula = "y ~ poly(x, 2)", se=F, size=1.5)+
    stat_smooth(data=filter(data, !(species_label %in% c("Colpidium striatum", "Dexiostoma campylum"))), aes(x=temperature, y=mean_log_major, colour=species_label, group=species_label, fill=species_label), method = "lm", formula = "y ~ x", se=F, size=1.5)+
       # stat_smooth(se=F, size=1.5, linetype="dashed")+
    xlab("temperature") + guides(colour="none", fill="none")  + scale_x_continuous(breaks=c(15,17,19,21,23,25)) + theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size=12), 
          axis.title = element_text(size=12)) + ggtitle("Size")
  p0b <- tag_facet(p0, x=Inf, hjust=1.5, open = c(""), close = c(")"), tag_pool = LETTERS) 
  
  # create composite figure on growth rate and max_bm
  p1 <- ggplot(data=subset(data), aes(x=temp_sc+20, y= estimate, colour=species_label, group=species_label, fill=species_label)) + 
    geom_point(size=1.5) + 
    facet_grid(species_label~., labeller = label_wrap_gen(width=10), scales = "free") +
    stat_smooth(data=filter(data, species_label %in% c("Colpidium striatum")), aes(x=temperature, y=estimate, colour=species_label, group=species_label, fill=species_label), method = "lm", formula = "y ~ poly(x, 2)", se=F, size=1.5)+
    stat_smooth(data=filter(data, !(species_label %in% c("Colpidium striatum"))), aes(x=temperature, y=estimate, colour=species_label, group=species_label, fill=species_label), method = "lm", formula = "y ~ x", se=F, size=1.5)+
#    stat_smooth(se=F, size=1.5, linetype="dashed")+
    ylab("Growth rate") + theme_bw() +  
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size=12), 
          axis.title = element_text(size=12), 
          strip.background.y = element_blank(), 
          strip.text.y = element_blank())  + xlab("temperature")+ guides(colour="none", fill="none")  + ggtitle("Demand")
  
  p2 <- ggplot(data=data, aes(x=temp_sc+20, y= AUC, colour=species_label, fill=species_label)) + geom_point(size=1.5) + 
    stat_smooth(data=filter(data, species_label %in% c("Colpidium striatum")), aes(x=temperature, y=AUC, colour=species_label, group=species_label, fill=species_label), method = "lm", formula = "y ~ poly(x, 2)", se=F, size=1.5)+
    stat_smooth(data=filter(data, !(species_label %in% c("Colpidium striatum"))), aes(x=temperature, y=AUC, colour=species_label, group=species_label, fill=species_label), method = "lm", formula = "y ~ x", se=F, size=1.5)+
    
    #stat_smooth(method = "lm", formula = "y ~ poly(x, 3)", se=F, linetype="dashed")+
    #stat_smooth(se=F, size=1.5, linetype="dashed")+
    facet_grid(species_label~., labeller = label_wrap_gen(width=10), scales="free") + ylab("Supply proxy") + xlab("temperature") + 
    guides(colour="none", fill="none")  + 
    theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size=12), 
          axis.title = element_text(size=12)) + ggtitle("Supply")
  
  
  p1 <- tag_facet(p1, x=Inf, hjust=1.5, open = c(""), close = c(")"), tag_pool = LETTERS[-c(1:6)]) 
  p2 <- tag_facet(p2, x=Inf, hjust=1.5, open = c(""), close = c(")"), tag_pool = LETTERS[-c(1:12)]) 
  
  
  p3 <- ggplot(data = data, aes(log_supply_over_log_demand, mean_log_major, group=species_label, colour=species_label)) +
    geom_point( size = 1.5) + 
    stat_smooth(method="lm", formula = 'y ~ x', se=F, fullrange = F)+
    #stat_smooth(se=F, size=1.5, linetype="dashed")+
     theme_bw() + facet_grid(species_label~., scales="free", labeller = label_wrap_gen(width=10)) +
    ylab("Mean log cell length") + 
    xlab("log(supply proxy) / log(demand)") + guides(colour="none", fill="none")  + 
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size=12), 
          axis.title = element_text(size=12)) + ggtitle("SD ratio")
  p3 <- tag_facet(p3, x=Inf, hjust=1.5, open = c(""), close = c(")"), tag_pool = LETTERS[-c(1:18)]) 
  p3 <- p3 + theme(strip.text = element_text(), strip.background = element_rect(fill = "grey", colour= "grey"))
  
  # p3 <- ggplot() +
  #   geom_point(data=data, aes(log_supply_over_log_demand, mean_log_major, group=species_label, colour=species_label), size = 1.5) +
  #   geom_line(data = RMA_results, aes(log_supply_over_log_demand, mean_log_major, group=species_label, colour=species_label), size=1.5) +
  #  geom_ribbon(data=RMA_results, aes(log_supply_over_log_demand, mean_log_major, ymin = conf_low, ymax = conf_high, fill=species_label, colour=NULL), alpha = 0.2) +
  #   geom_text(data=RMA_results, mapping = aes(-Inf, Inf,hjust = -0.1, vjust = 2, label = lab), colour="black", parse=F, check_overlap = TRUE) +
  #   theme_bw() + facet_grid(species_label~., scales="free", labeller = label_wrap_gen(width=10)) +
  #   ylab("Mean log cell length") + 
  #   xlab("log(supply proxy) / log(demand)") + guides(colour=F, fill=F) + 
  #   theme(axis.text = element_text(size=12), 
  #         axis.title = element_text(size=12))
  # p3 <- tag_facet(p3, x=Inf, hjust=1.5, open = c(""), close = c(")"), tag_pool = LETTERS[-c(1:18)]) 
  # p3 <- p3 + theme(strip.text = element_text(), strip.background = element_rect(fill = "grey", colour= "grey"))
  
  
  cowplot::plot_grid(p0b, p1, p2, p3, ncol=4)
  ggsave(here::here("output/figures/main/Figure_1.jpg"), width = 12, height=8)
  
}









create_poly_SD_data <- function(params, data, supply_proxy = "AUC"){

  #browser()
  
  pop_dyn_poly <- data %>% filter(day <= 40) %>% group_by(predicted_species, combination, richness, temperature, day, replicate, microcosmID) %>%
    dplyr::summarize(response_var = sum(bvol, na.rm = T))
  
  dd_auc_poly <- pop_dyn_poly %>% arrange(combination, richness, predicted_species, temperature, replicate, day) %>% 
    group_by(microcosmID, combination, richness, predicted_species, temperature, replicate) %>% nest()
  
  dd_auc_poly$AUC <- unlist(lapply(1:nrow(dd_auc_poly), function(x) area_under_curve(response_var = dd_auc_poly$data[[x]]$response_var, time=dd_auc_poly$data[[x]]$day)))
  
  dd_auc_poly$max_bm <- unlist(lapply(1:nrow(dd_auc_poly), function(x) max( dd_auc_poly$data[[x]]$response_var, na.rm=T)))
  
  dd_auc_poly_size <- data %>% filter(day <= 40) %>%  group_by(combination, richness, predicted_species, temperature, replicate) %>% dplyr::summarize(mean_log_major = mean(log(Major), na.rm=T))
  
  dd_auc_poly_size <- merge(dd_auc_poly, dd_auc_poly_size)
  
  mono_growth_rates <- subset(params, term == "r") %>% group_by(predicted_species, temperature) %>% dplyr::summarize(estimate = mean(estimate, na.rm=T))
  
  # merge with the growth rates measured in monocultures
  dd_auc_poly_size <- merge(dd_auc_poly_size, mono_growth_rates)
  dd_auc_poly_size <- dd_auc_poly_size %>% mutate(log_supply_over_log_demand = log(get(supply_proxy)) - log(estimate))
  dd_auc_poly_size$sp_rich <- paste0(dd_auc_poly_size$predicted_species, "_", dd_auc_poly_size$richness)
  
  # add full names for species
  dd_auc_poly_size <- dd_auc_poly_size %>% mutate(species_label =
                                                    case_when(predicted_species == "Colp" ~ "Colpidium striatum",
                                                              predicted_species == "Dexio" ~ "Dexiostoma campylum",
                                                              predicted_species == "Loxo" ~ "Loxocephalus sp.",
                                                              predicted_species == "Para" ~ "Paramecium caudatum",
                                                              predicted_species == "Spiro" ~ "Spirostomum teres",
                                                              predicted_species == "Tetra" ~ "Tetrahymena thermophila"),
                                                  incubator = case_when(
                                                    as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
                                                    as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
                                                    as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
                                                    as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
                                                    as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
                                                    as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
                                                    as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
                                                    as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
                                                    as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
                                                    as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
                                                    as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
                                                    as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
                                                  ))
  
  return(dd_auc_poly_size)
  
  
}




create_poly_output_figure <- function(data){

  data <- data %>% mutate(temperature = case_when(
    round(temp_sc, 2) == round(-5.0024073, 2) ~ 15,
    round(temp_sc, 2) == round(-3.0024073, 2) ~ 17,
    round(temp_sc, 2) == round(-1.0024073, 2) ~ 19,
    round(temp_sc, 2) == round(0.9975927, 2) ~ 21,
    round(temp_sc, 2) == round(2.9975927, 2) ~ 23,
    round(temp_sc, 2) == round(4.9975927, 2) ~ 25
  ))

    polyplot <- ggplot(data=subset(data, richness_fac != 1)) + 
    geom_point(aes(x=temperature, y=log(mean_size), colour=as.numeric(richness_fac), group=interaction(species_label, combination))) +  
    geom_line(aes(x=temperature, y=(pred_size), group=interaction(species_label, combination), colour=as.numeric(richness_fac)))  + theme_bw() +     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=12)) +
    facet_wrap(~species_label, scale = "free", ncol = 6, labeller = label_wrap_gen(width=10)) +  scale_colour_viridis_c() +
    ylab("Mean log cell length") +
    xlab("Temperature") + guides(colour = guide_colourbar(title="richness")) + scale_x_continuous(breaks=c(15,20,25))
  
  polyplot2 <- tag_facet(polyplot, x=Inf, hjust=1.5, open = c(""), close = c(")"), tag_pool = 1:6) 
  polyplot2 <- polyplot2 + theme(strip.text = element_text(), strip.background = element_rect(fill = "grey", colour= "grey"))
  
  polyplot2
  #browser()
  #cowplot::ggsave2(here::here("output/figures/main/Figure_3.jpg"), width = 12, height=4)
}



model_plot_SD_poly <- function(data, figure){
  
#data <- bind_rows(poly_SD_data)
data <- data %>% group_by(predicted_species) %>% 
  mutate(log_supply_over_log_demand_c = log_supply_over_log_demand-mean(log_supply_over_log_demand),
         rich_sc = richness - mean(richness, na.rm=T))
data$rich_fac <- as.factor(data$richness)
#dd <- bind_rows(data)
#dd$species <- rep(c("Colp", "Dexio", "Loxo", "Para", "Spiro", "Tetra"), each = 27)
#i <- "Dexio"

counter <- 1
models <- vector("list", 6)
models[[counter]] <- list()


for (i in c("Colp", "Dexio", "Loxo", "Para", "Spiro", "Tetra")){

selected_species <- i

#pdf(here::here(paste0("output/figures/not_for_paper/Figure_SX_", i,".pdf")), width = 10, height=8)

#data$comb_rich <- paste0(data$combination,"_", data$richness)

 sp_data <- subset(data, predicted_species == selected_species & richness > 1) #%>% mutate(center_log_supply_over_log_demand = log_supply_over_log_demand - mean(log_supply_over_log_demand, na.rm=T))
 
 mod1 <- lmer(mean_log_major ~ log_supply_over_log_demand_c * rich_sc + (1|incubator) + (1|combination), 
              data=sp_data, REML=T)
 
# mod1 <- lmer(mean_log_major ~ log_supply_over_log_demand + (1 | rich_fac:combination), 
#              data=subset(data, predicted_species == selected_species & richness > 1), REML=T)



# print(dwplot(x=list(mod1))+geom_vline(xintercept=0,lty=2) + ggtitle("Mixed model (random intercepts)"))
# print(plot(mod1,col=as.factor(sp_data$combination),id=0.05))
#   # cols2 <- c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
#   # cols3 <- c("#edf8e9", "#bae4b3",   "#74c476", "#31a354", "#006d2c")
#   # cols4 <- c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#3f007d")
#   # cols5 <- c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
#   # cols6 <- c("#252525")
#   # cols <- c(cols2, cols3, cols4, cols5, cols6)
#   
#   tt <- tidy(mod1,effects="ran_vals")
#   #  tt$richness <- str_count(tt$level, ',') +1
#   #  combis <- tt %>% arrange(richness) %>% select(level)
#   tt$level2 <- factor(tt$level, levels=tt$level[order(tt$estimate)], ordered=TRUE)
#   
#   #  cols <- setNames(cols, combis$level)  
#   
#   print(ggplot(tt,aes(level2,estimate))+
#           geom_pointrange(aes(ymin=estimate-1.96*std.error,
#                               ymax=estimate+1.96*std.error))+
#           theme(legend.position="bottom")+
#           coord_flip() + ggtitle("Random effects (intercept)")) 
  
  
  # mod2 <- lmer(mean_log_major ~ log_supply_over_log_demand * poly(richness, 1)  + (0 + log_supply_over_log_demand | combination), 
  #              data=sp_data, REML=T)
  # summary(mod2)
  # 
  # 
  # print(dwplot(mod2)+geom_vline(xintercept=0,lty=2) + ggtitle("Mixed model (random slopes)"), control = lmerControl(
  #   optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  # print(plot(mod2,col=as.factor(sp_data$combination),id=0.05))
  # 
  # 
  # tt <- tidy(mod2,effects="ran_vals")
  # tt$level2 <- factor(tt$level, levels=tt$level[order(tt$estimate)], ordered=TRUE)
  # 
  # print(ggplot(tt,aes(level2,estimate))+
  #   geom_pointrange(aes(ymin=estimate-1.96*std.error,
  #                       ymax=estimate+1.96*std.error))+
  #     theme(legend.position="bottom")+
  #     coord_flip() + ggtitle("Random effects (intercept)")) 
  
  # mod3 <- lmer(mean_log_major ~ log_supply_over_log_demand * poly(richness, 2) + (1 + log_supply_over_log_demand | combination), 
  #              data=subset(data, predicted_species == selected_species & richness > 1), REML=F, control = lmerControl(
  #                optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  
  # mod3 <- lmer(mean_log_major ~ log_supply_over_log_demand * poly(richness, 1) + (1 + log_supply_over_log_demand | combination), data=sp_data, REML=T, control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  # #mod3 <- mod2
  # summary(mod3)
  # 
  # ggplot(data=sp_data, aes(x=mean_log_major)) + geom_histogram(binwidth = 0.0001)
  # 
  # #ss <- getME(mod3,c("theta","fixef"))
  # #mod4 <- update(mod3,start=ss,control=lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  # 
  # # check singularity
  # #tt <- getME(mod1,"theta")
  # #ll <- getME(mod1,"lower")
  # #min(tt[ll==0])
  # 
  # print(dwplot(mod3)+geom_vline(xintercept=0,lty=2) + ggtitle("Mixed model (random regression)"))
  # print(plot(mod3,col=as.factor(sp_data$combination),id=0.05))
  
  
  #mod1 <- lm(mean_log_major ~ log_supply_over_log_demand * combination + as.factor(richness), data=data)  
  # browser()
  #  library(rstanarm)
  #  
  #  blmer <- stan_lmer(mean_log_major ~ log_supply_over_log_demand * rich_sc + (1 + log_supply_over_log_demand | combination), 
  #                     data=data, iter=10000)
  #  summary(blmer, pars = "beta", digits = 3)
  #  pairs(blmer)
  
 # print(AIC(mod1, mod2, mod3))
  
    models[[counter]] [['Intercept']] <- mod1
  #  models[[counter]] [['Slope']] <- mod2
  #  models[[counter]] [['Intercept_Slope']] <- mod3
    
  dataf <- sp_data
  dataf$pred_int <- predict(mod1, re.form=~(1|combination)) 
  #dataf$pred_slope <- predict(mod2) 

  # print(ggplot(data=dataf, aes(y=mean_log_major, x=log_supply_over_log_demand, group=combination, colour=richness)) + 
  #         geom_line(aes(y=pred_slope,x=log_supply_over_log_demand))+
  #         geom_line(aes(y=pred_int,x=log_supply_over_log_demand), linetype="dashed")+
  #         geom_point() + facet_wrap(~richness, shrink = T) + guides(colour=F)) + ggtitle(selected_species)
  # 
  # print(ggplot(data=dataf, aes(y=mean_log_major, x=log_supply_over_log_demand, group=combination, colour=combination)) + 
  #   geom_line(aes(y=pred_slope,x=log_supply_over_log_demand))+
  #   geom_line(aes(y=pred_int,x=log_supply_over_log_demand), linetype="dashed")+
  #   geom_point() + facet_wrap(~combination, shrink = T) + guides(colour=F)) + ggtitle(selected_species)

  plt <- ggplot(data=dataf, aes(y=mean_log_major, x=log_supply_over_log_demand_c, group=combination, colour=richness)) + 
          geom_line(aes(y=pred_int,x=log_supply_over_log_demand))+ scale_colour_viridis_c() +  theme_bw() +
    ylab("Mean log cell length") + 
    xlab(expression(over("log(supply proxy)", "log(growth rate)"))) 
  
  assign(paste0("plt_", selected_species), plt)
  
#dev.off()
counter <- counter+1
}

# give names to models
names(models) <- c("Colpidium striatum", "Dexiostoma campylum", "Loxocephalus sp.", "Paramecium caudatum", "Spirostomum teres", "Tetrahymena thermophila")
#print(models)

cm <- c('(Intercept)' = 'Constant',
        'log_supply_over_log_demand' = "SD ratio",
        'rich_sc' = 'Scaled richness',
        'log_supply_over_log_demand:rich_sc' = 'Temperature : richness',
        'cor__(Intercept).temp_sc' = "Intercept-slope correlation",
        'sd__(Intercept)' = 'SD Intercept',
        'sd__Observation' = "SD observation",
        'sd__temp_sc' = "SD temperature slope"
)

options(modelsummary_get = "easystats")
# msummary(lapply(models, "[[", 1), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC',
#          output = here::here("output/tables/mixed_model_table_poly_SD_random_int.docx"))

# msummary(lapply(models, "[[", 2), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC',
#          output = here::here("output/tables/mixed_model_table_poly_SD_random_slope.docx"))
# 
# msummary(lapply(models, "[[", 3), statistic = 'conf.int', conf_level = .95, fmt = "%.4f", gof_omit = 'DF|Deviance|Log.Lik.|REMLcrit|BIC',
#          output = here::here("output/tables/mixed_model_table_poly_SD_random_intercept_slope.docx"))

print(figure / (plt_Colp + ggtitle(c("Colpidium striatum")) + xlab("") + 
                                        plt_Dexio + ggtitle(c("Dexiostoma campylum")) + xlab("") + ylab("") + 
                                        plt_Loxo + ggtitle(c("Loxocephalus sp.")) + xlab("") + ylab("") + 
                                        plt_Para + ggtitle(c("Paramecium caudatum")) + 
                                        plt_Spiro + ggtitle(c("Spirostomum teres")) + ylab("") + 
                                        plt_Tetra+ ylab("") + ggtitle(c("Tetrahymena thermophila")) +
                                       plot_layout(ncol = 3, guides = 'auto', axis_titles = 'collect_x')) + 
             plot_layout(guides = 'collect', axis_titles = 'collect_y') + plot_annotation(tag_levels = 'A'))
ggsave("output/figures/main/Figure_4.jpg", height=10,width=10)

}



