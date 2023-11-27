#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status.

# Function to log messages
log() {
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1"
}

# Function to clean up in case of an error
cleanup() {
  log "An error occurred. Cleaning up..."
  rm -f /path/to/a/temporary/file.txt
  exit 1
}

# Function to send a notification (placeholder function)
send_notification() {
  log "Sending notification about failure..."
  # Add actual notification code here (e.g., send email)
}

# Trap signals and errors
trap cleanup ERR SIGINT SIGTERM

# Main loop
cat /path/to/list/of/reads.txt | while read lines; do
  touch /path/to/a/temporary/file.txt
  echo "$lines" >> /path/to/a/temporary/file.txt

  cat /path/to/a/temporary/file.txt | while read line; do
    # Attempt to run snakemake up to 3 times if it fails
    local max_retries=3
    for ((i=1; i<=max_retries; i++)); do
      log "Starting snakemake for $line (Attempt $i of $max_retries)"
      if time snakemake --cores 16 -F; then
        log "Snakemake completed successfully for $line"
        break
      else
        log "Snakemake failed for $line (Attempt $i of $max_retries)"
        if [ "$i" -eq "$max_retries" ]; then
          log "Max retries reached for $line. Moving to the next individual."
          send_notification
          break
        fi
        sleep 10 # Wait for 10 seconds before retrying
      fi
    done

    echo "$line" >> /path/to/list/of/analyzed/individuals.txt
  done

  rm /path/to/a/temporary/file.txt
done

log "Pipeline completed."