library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(pROC)
library(glmnet)
library(ggplot2)
library(hoopR)

################### Settings
min_minutes <-8 
min_games <- 20
nclus <- 35
current_year <- 2022

################### Read the data
nba_pbp <- hoopR::load_nba_player_box(2017:2022)
box_df <- 
  filter(nba_pbp, as.numeric(season_type)>1 & min != "--") %>%
  separate(fg3, into=c("fg3_made", "fg3_attempts"), sep="-") %>%
  separate(fg, into=c("fg_made", "fg_attempts"), sep="-") %>%
  separate(ft, into=c("ft_made", "ft_attempts"), sep="-") %>%
  mutate(row=row_number(),
         season=as.numeric(season),
         season_type=as.numeric(season_type),
         game_date=as.Date(game_date)-1,
         game_id=as.numeric(game_id),
         team_id=as.numeric(team_id),
         athlete_id=as.numeric(athlete_id),
         min=as.numeric(min), 
         played_min=if_else(min>min_minutes, 1, 0),
         s=min/48,
         oreb=as.numeric(oreb)/s, 
         dreb=as.numeric(dreb)/s, 
         stl=as.numeric(stl)/s, 
         ast=as.numeric(ast)/s,
         blk=as.numeric(blk)/s, 
         to=as.numeric(to)/s,
         pts=as.numeric(pts), 
         ft_made=as.numeric(ft_made), 
         fg_made=as.numeric(fg_made), 
         fg3_made=as.numeric(fg3_made), 
         ft_attempts=as.numeric(ft_attempts), 
         fg_attempts=as.numeric(fg_attempts), 
         fg3_attempts=as.numeric(fg3_attempts),
         attempts=(ft_attempts+fg_attempts+fg3_attempts)/s,
         plus_minus_n=if_else(grepl("+", plus_minus), as.numeric(gsub("+", "", plus_minus)), -as.numeric(gsub("-", "", plus_minus)))) %>%
  filter(!(team_abbreviation %in% c("LEB", "USA", "WORLD", "GIA", "DUR"))) %>%
  group_by(game_id) %>%
  mutate(home=if_else(row==max(row), 1, 0)) %>%
  group_by(game_id, team_id) %>%
  mutate(home=max(home),
         team_pts=sum(pts), 
         share_min=min/sum(min), 
         current_playoffs=if_else(season==current_year & season_type==3, 1, 0)) %>%
  group_by(athlete_id) %>%
  mutate(games_played_ever=sum(played_min)) %>%
  arrange(game_date, game_id, team_id) %>%
  dplyr::select(home, games_played_ever, current_playoffs, season, season_type, athlete_id, athlete_display_name, team_id, team_abbreviation, game_id, game_date, team_pts, oreb, attempts, dreb, stl, ast, min, blk, to, pts, ft_made, ft_attempts, fg_made, fg_attempts, fg3_made, fg3_attempts, share_min, plus_minus_n) %>%
  ungroup() %>%
  filter(!between(game_date, as.Date("2020-07-09"), as.Date("2020-10-11")))
  


################### Generate clusters
### First get the normalized stats
players_df_collapsed <- 
  filter(box_df, games_played_ever>min_games & min>min_minutes) %>%
  group_by(athlete_id) %>%
  summarise(across(oreb:plus_minus_n, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(fg3_percent=ifelse(fg3_attempts>0, fg3_made/fg3_attempts, 0),
         fg_percent=ifelse(fg_attempts>0, fg_made/fg_attempts, 0), 
         ft_percent=ifelse(ft_attempts>0, ft_made/ft_attempts, 0)) %>%
  dplyr::select(-min, -fg3_made, -ft_made,  -ft_attempts, -fg_made, -pts, -fg3_attempts, -fg_attempts, -athlete_id) %>%
  ungroup() %>%
  mutate(across(everything(), ~ scale(.)[,1])) 

set.seed(2022)
km <- kmeans(players_df_collapsed, centers=nclus, nstart=50, iter.max=100)
centroids <- km$centers

### Check the choice of K
data.frame(k=rep(1:50), 
           wss=sapply(1:50, function(k){kmeans(players_df_collapsed, k, nstart=50, iter.max = 15)$tot.withinss/km$totss})) %>%
  ggplot(aes(x=k, y=wss)) + geom_line()
  
### Function to map clusters to athletes
closest.cluster <- function(x) {
  cluster.dist <- apply(centroids, 1, function(y) sqrt(sum((x-y)^2)))
  return(which.min(cluster.dist)[1])
}

### Normalized stats by season so athletes can change clusters across seasons
players_df_by_season <- 
  group_by(box_df, athlete_id, athlete_display_name, season) %>%
  filter(games_played_ever>min_games & min>min_minutes) %>%
  summarise(across(oreb:plus_minus_n, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(attempts=ft_attempts+fg_attempts+fg3_attempts,
         fg3_percent=ifelse(fg3_attempts>0, fg3_made/fg3_attempts, 0),
         fg_percent=ifelse(fg_attempts>0, fg_made/fg_attempts, 0), 
         ft_percent=ifelse(ft_attempts>0, ft_made/ft_attempts, 0)) %>%
  ungroup()

players_df_by_season_std <- 
  group_by(players_df_by_season, season) %>%
  dplyr::select(-min, -fg3_made, -ft_made,  -ft_attempts, -fg_made, -pts, -fg3_attempts, -fg_attempts, -athlete_id, -athlete_display_name) %>%
  mutate(across(everything(), ~ scale(.)[,1])) %>%
  ungroup() %>%
  dplyr::select(-season)

### Assign the clusters
clusters <- data.frame(cbind(apply(players_df_by_season_std, 1, closest.cluster), players_df_by_season$athlete_id, players_df_by_season$season, players_df_by_season$athlete_display_name), stringsAsFactors = FALSE)
names(clusters) <- c("cluster", "athlete_id", "season", "athlete_display_name")  
clusters$season <- as.numeric(clusters$season)
clusters$cluster <- as.numeric(clusters$cluster)
clusters$athlete_id <- as.numeric(clusters$athlete_id)


## Check the  counts
dist <- group_by(clusters, cluster) %>%
  summarise(n=n())

set.seed(2022)
################### Set up data for modeling by calculating the net differences in cluster distributions
box_df_clusters <- left_join(box_df, dplyr::select(clusters, -athlete_display_name), by=c("athlete_id", "season")) %>%
  replace_na(list(cluster=0)) %>%
  group_by(current_playoffs, season, season_type, game_id, game_date, team_abbreviation, team_id, team_pts, cluster, home) %>%
  summarise(min=sum(min)) %>%
  arrange(game_id, team_pts) %>%
  group_by(game_id) %>%
  mutate(win=if_else(dense_rank(team_pts)==2, 1, 0)) %>%
  ungroup() %>%
  pivot_wider(id_cols=c("game_id", "team_id", "game_date", "team_abbreviation", "win", "team_pts", "season", "season_type", "current_playoffs", "home"), names_from = cluster, values_from = min, names_prefix = "min_", values_fill = 0) %>%
  group_by(game_id) %>%
  mutate(r = sample(100,1), 
         team_indicator=case_when(r<50 & row_number()==1 ~ 2, 
                                  r<50 & row_number()==2 ~ 1,
                                  TRUE ~ as.numeric(row_number()))) %>%
  ungroup() %>%
  dplyr::select(-r) %>%
  arrange(game_id, team_indicator) %>%
  mutate(total_min=rowSums(across(starts_with("min"))), 
         m=if_else(team_indicator==1, 1, -1), 
         label=if_else(team_indicator==1 & win==1, 1, 0), 
         home_team=if_else(team_indicator==1 & home==1, 1, 0)) %>% 
  mutate_at(vars(starts_with("min")), ~(. / total_min*m)) %>%
  arrange(game_date, game_id, team_id) %>%
  ungroup()

box_df_clusters_agg <-
  group_by(box_df_clusters, game_id, season, season_type, current_playoffs) %>%
  summarise_at(vars(starts_with("min"), label, home_team), sum) %>%
  ungroup()

################### Fit the model (excluding the modified 2022 playoff data for model estimation)
t <- filter(box_df_clusters_agg)
train <- sample_frac(t, .7) 
valid <- filter(t, !(game_id %in% train$game_id))
rm(t)

Y <- train$label
Y_VALID <- valid$label
X <- model.matrix(as.formula(Y ~ .), dplyr::select(train, starts_with("min"), home_team))
X_VALID <- model.matrix(as.formula(Y_VALID ~ .), dplyr::select(valid, starts_with("min"), home_team))


model <- cv.glmnet(y=Y, x=X, family="binomial", parallel=FALSE, nfolds=10)
c <- as.matrix(coef(model, s=model$lambda.1se))

### Check the AUC
p_valid <- as.numeric(1/(1+exp(-X_VALID%*%c[-1])))
auc(roc(Y_VALID, p_valid))

p <- as.numeric(1/(1+exp(-X%*%c[-1])))
auc(roc(Y, p))

################### Run the playoffs
### Get the most recent cluster distributions from the modified playoff data
playoffs <- 
  left_join(box_df, dplyr::select(clusters, -athlete_display_name), by=c("athlete_id", "season")) %>%
  replace_na(list(cluster=0)) %>%
  filter(athlete_display_name != "Joel Embiid") %>% # injured
  filter(athlete_display_name != "Khris Middleton") %>% # injured
  group_by(game_id, game_date, team_abbreviation, team_id, cluster, current_playoffs) %>%
  summarise(min=sum(min)) %>%
  ungroup() %>%
  pivot_wider(id_cols=c("game_id", "team_id", "game_date", "team_abbreviation", "current_playoffs"), names_from = cluster, values_from = min, names_prefix = "min_", values_fill = 0) %>%
  mutate(total_min=rowSums(across(starts_with("min")))) %>% 
  mutate_at(vars(starts_with("min")), ~(. / total_min)) %>%
  arrange(game_date, game_id, team_id) %>%
  filter(current_playoffs==1) %>%
  dplyr::select(team_abbreviation, starts_with("min")) %>%
  ungroup()

### Helper function to get the probability from the elastic net model and simulate the 7 game series
matchup <- function(team1, team2, played){
  team1_stats <- filter(playoffs, team_abbreviation==team1)  %>%
    summarise(across(starts_with("min"), ~ mean(abs(.x), na.rm = TRUE))) %>%
    mutate(total_min=rowSums(across(starts_with("min")))) %>%
    mutate_at(vars(starts_with("min")), ~(. / total_min)) %>%
    mutate(total_min=rowSums(across(starts_with("min"))))

  #print(team1_stats$total_min)
  
  team2_stats <- filter(playoffs, team_abbreviation==team2)  %>%
    summarise(across(starts_with("min"), ~ mean(abs(.x), na.rm = TRUE))) %>%
    mutate(total_min=rowSums(across(starts_with("min")))) %>%
    mutate_at(vars(starts_with("min")), ~(. / -total_min)) %>%
    mutate(total_min=rowSums(across(starts_with("min"))))
  
  #print(team2_stats$total_min)

  input1 <- bind_rows(team1_stats, team2_stats) %>%
    summarise_at(vars(starts_with("min")), sum) %>%
    mutate(home_team=1) %>%
    ungroup()

  input2 <- 
    mutate(input1, home_team=0) 
  
  Y <- 0
  X1 <- model.matrix(as.formula(Y ~ .), dplyr::select(input1, starts_with("min"), home_team))
  p1 <- as.numeric(1/(1+exp(-X1[,row.names(c)[-1]]%*%c[-1])))
  X2 <- model.matrix(as.formula(Y ~ .), dplyr::select(input2, starts_with("min"), home_team))
  p2 <- as.numeric(1/(1+exp(-X2[,row.names(c)[-1]]%*%c[-1])))
  p <- c(p1, p1, p2, p2, p1, p2, p1)
  
  n1 <- 0
  n2 <- 0
  for (k in 1:25000){
    nn1 <- 0
    nn2 <- 0
    for (g in 1:7){
      if (played[g]==1){
        binomial <- 1
      } else if (played[g]==-1){
        binomial <- 0
      } else{
        binomial <- as.numeric(rbinom(n=1, size=1, prob=p[g]))  
      }
      if (binomial==1){nn1 <- nn1+1}
      else {nn2 <- nn2+1}
    }
    if (nn1>nn2){
      n1 <- n1+1
    } else{
      n2 <- n2+1
    }
  }
  prob_win <- n1/(n1+n2)
  return(prob_win)
}

contribution <- function(team1, team2){
  team1_stats <- filter(playoffs, team_abbreviation==team1)  %>%
    summarise(across(starts_with("min"), ~ mean(abs(.x), na.rm = TRUE))) %>%
    mutate(total_min=rowSums(across(starts_with("min")))) %>%
    mutate_at(vars(starts_with("min")), ~(. / total_min)) %>%
    mutate(total_min=rowSums(across(starts_with("min"))))
  
  team2_stats <- filter(playoffs, team_abbreviation==team2)  %>%
    summarise(across(starts_with("min"), ~ mean(abs(.x), na.rm = TRUE))) %>%
    mutate(total_min=rowSums(across(starts_with("min")))) %>%
    mutate_at(vars(starts_with("min")), ~(. / -total_min)) %>%
    mutate(total_min=rowSums(across(starts_with("min"))))
  
  inputx <- bind_rows(team1_stats, team2_stats) %>%
    summarise_at(vars(starts_with("min")), sum) %>%
    mutate(home_team=0) %>%
    ungroup()

  Y <- 0
  XX <- model.matrix(as.formula(Y ~ .), dplyr::select(inputx, starts_with("min"), home_team))
  print(XX[,row.names(c)[-1]])
  d <- data.frame(row.names(c)[-1], XX[,row.names(c)[-1]]*c[-1])
  names(d) <- c("variable", "contribution") 
  return(d)

}

View(contribution("GS","BOS"))

### Run the games for round 2

matchup("GS","BOS", c(-1,1,-1,1,1,0,0))

