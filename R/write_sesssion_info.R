write_session_info <- function(dir){
    writeLines(utils::capture.output(sessioninfo::session_info()),
               paste(dir, "sessionInfo.txt", sep = "/"))
    }
