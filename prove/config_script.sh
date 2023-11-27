#!/bin/bash

echo "Welcome to the RePaRe configuration setup."

# Function to ask for a path, validate it, and write to config
ask_path() {
    local prompt=$1
    local config_key=$2
    local path
    read -p "$prompt" path

    # Validate the path
    if [ -z "$path" ]; then
        echo "You must specify a path."
        ask_path "$prompt" "$config_key"
    else
        echo "Setting $config_key to $path"
    # Write to config
    else
        echo "Setting $config_key to $path"
        echo "$config_key=$path" >> "$HOME/.repare.cfg"
}

# Ask for paths
ask_path "Enter the absolute path to your list of samples: " "list_path_first"
ask_path "Enter the absolute path to your working directory: " "working_dir"
# ... other configurations

echo "Configuration setup is complete."