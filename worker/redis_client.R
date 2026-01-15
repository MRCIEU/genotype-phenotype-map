library(futile.logger)

process_gwas <- 'process_gwas'
process_gwas_dlq <- glue::glue('{process_gwas}_dlq')
delete_gwas_queue <- 'delete_gwas'

connect_to_redis <- function(max_retries = 5, initial_delay = 2) {
  host <- Sys.getenv("REDIS_HOST", "redis")
  port <- as.numeric(Sys.getenv("REDIS_PORT", 6379))
  
  for (attempt in 1:max_retries) {
    tryCatch({
      conn <- redux::hiredis(host = host, port = port)
      conn$PING()
      flog.info(paste("Successfully connected to Redis on attempt", attempt))
      return(conn)
    }, error = function(e) {
      if (attempt < max_retries) {
        delay <- initial_delay * (2 ^ (attempt - 1))
        flog.warn(paste("Failed to connect to Redis on attempt", attempt, "of", max_retries, "- retrying in", delay, "seconds. Error:", e$message))
        Sys.sleep(delay)
      } else {
        flog.error(paste("Failed to connect to Redis after", max_retries, "attempts. Error:", e$message))
        stop(paste("Could not connect to Redis after", max_retries, "attempts"))
      }
    })
  }
}

send_to_dlq <- function(redis_conn, message) {
  redis_conn$LPUSH(process_gwas_dlq, message)
}

get_from_process_queue <- function(redis_conn) {
  redis_conn$BRPOP(process_gwas, timeout = 0.1)
}

get_from_delete_queue <- function(redis_conn) {
  redis_conn$BRPOP(delete_gwas_queue, timeout = 0.1)
}