#! /usr/bin/bash -x

echo "Creating a ABCSMC skeleton. WARNING: may overwrite existing files!"

# This script will create a skeleton for a new abc-project
# Create the target directory for the skeleton
mkdir -p $1

cp -r template/* $1