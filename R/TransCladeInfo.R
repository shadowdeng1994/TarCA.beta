#' Transform clade infomation into a dataframe
#'
#' @param CCChar A character of clade infomation.
#'
#' @return A dataframe with the columns of clade size and frequency.
#' @export
#'
#' @examples
#' library(TarCA)
#' tmp.cladesize <- "1 (12), 2 (1), 3 (2)"
#' fun.TransCladeInfo(tmp.cladesize)
#'

fun.TransCladeInfo <- function(CCChar){
  stringr::str_split(CCChar,",") %>% unlist %>%
    tbl_df %>%
    mutate(value=value %>% gsub("^ +","",.) %>% gsub("[()]","",.)) %>%
    separate(value,c("Size","Freq")," ") %>%
    mutate(Size=Size %>% as.numeric,Freq=Freq %>% as.numeric)
}
