#!/usr/bin/env bash

# Define configuration file names
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG_ENV_FNAME="config_env.sh"
UNCONFIG_ENV_FNAME="unconfig_env.sh"
CONFIG_ENV_PATH="$SCRIPT_DIR/$CONFIG_ENV_FNAME"
UNCONFIG_ENV_PATH="$SCRIPT_DIR/$UNCONFIG_ENV_FNAME"

# Determine the project root directory
ENV_FILE="./.env"
if [ ! -f "$ENV_FILE" ]; then
  echo ".env file not found in $SCRIPT_DIR"
  exit 1
fi

# Load the .env file
set -a
. "$ENV_FILE"
set +a

# Resolve LRGSG_ROOT to the directory where the .env file is located
LRGSG_ROOT="$(cd "$(dirname "$ENV_FILE")" && pwd)"

# Export LRGSG_ROOT explicitly
export LRGSG_ROOT

# Update paths to use the inferred LRGSG_ROOT
for path in $(grep -o '^[^#]*' "$ENV_FILE" | grep '=' | cut -d '=' -f 1); do
    if [ "$path" != "LRGSG_ROOT" ]; then
        eval "export $path"
    fi
done

# Generate config_env.sh
cat << EOF > $CONFIG_ENV_PATH
#!/bin/sh

# Export directories
export LRGSG_ROOT="$LRGSG_ROOT"
EOF

# for path in $(grep -o '^[^#]*' "$ENV_FILE" | grep '=' | cut -d '=' -f 1); do
#     value=$(eval echo \$$path)
#     echo "export $path=\"$value\"" >> $CONFIG_ENV_PATH
# done

for path in $(grep -o '^[^#]*' "$ENV_FILE" | grep '=' | cut -d '=' -f 1); do
    if [ "$path" != "LRGSG_ROOT" ]; then
        raw=$(grep "^$path=" "$ENV_FILE" | cut -d '=' -f2-)
        # replace any literal ${LRGSG_ROOT} or $LRGSG_ROOT with the computed one
        value="${raw//\$\{LRGSG_ROOT\}/$LRGSG_ROOT}"
        value="${value//\$LRGSG_ROOT/$LRGSG_ROOT}"
        echo "export $path=\"$value\"" >> $CONFIG_ENV_PATH
    fi
done


# Generate unconfig_env.sh
cat << EOF > $UNCONFIG_ENV_PATH
#!/bin/sh

# Unset directories
EOF

for path in $(grep -o '^[^#]*' "$ENV_FILE" | grep '=' | cut -d '=' -f 1); do
    echo "unset $path" >> $UNCONFIG_ENV_PATH
done

echo "Scripts $CONFIG_ENV_FNAME and $UNCONFIG_ENV_FNAME have been generated."



# #!/bin/sh

# # Define configuration file names
# SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# CONFIG_ENV_FNAME="config_env.sh"
# UNCONFIG_ENV_FNAME="unconfig_env.sh"
# CONFIG_ENV_PATH="$SCRIPT_DIR/$CONFIG_ENV_FNAME"
# UNCONFIG_ENV_PATH="$SCRIPT_DIR/$UNCONFIG_ENV_FNAME"
# # Source the paths.sh file
# . "$SCRIPT_DIR/paths.sh"

# # Generate config_env.sh
# cat << EOF > $CONFIG_ENV_PATH
# #!/bin/sh

# # Determine the project root directory based on the location of this script
# LRGSG_ROOT="\$(cd "\$(dirname "\$0")" && pwd)"

# # Export directories
# export LRGSG_ROOT
# EOF

# for path in "${paths[@]}"; do
#     echo "export ${path}" >> $CONFIG_ENV_PATH
# done

# # Generate unconfig_env.sh
# cat << EOF > $UNCONFIG_ENV_PATH
# #!/bin/sh

# # Unset directories
# EOF

# for path in "${paths[@]}"; do
#     var_name=$(echo "${path}" | cut -d '=' -f 1)
#     echo "unset ${var_name}" >> $UNCONFIG_ENV_PATH
# done

# echo "Scripts $CONFIG_ENV_FNAME and $UNCONFIG_ENV_FNAME have been generated."
