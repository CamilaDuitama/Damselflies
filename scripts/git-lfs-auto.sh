#!/bin/bash

# Set the size threshold to 100MB (in bytes)
SIZE_THRESHOLD=$((100 * 1024 * 1024))

# Find all files larger than the threshold
find . -type f -size +${SIZE_THRESHOLD}c | while read file; do
  # Skip .git directory
  if [[ $file == ./.git/* ]]; then
    continue
  fi
  
  # Get the relative path
  relative_path=${file#./}
  
  # Add to Git LFS
  git lfs track "$relative_path"
  
  # Add to .gitignore
  echo "$relative_path" >> .gitignore
  
  echo "Added $relative_path to Git LFS and .gitignore"
done

# Remove duplicate lines from .gitignore
sort -u .gitignore > .gitignore.tmp && mv .gitignore.tmp .gitignore

# Remove duplicate lines from .gitattributes
sort -u .gitattributes > .gitattributes.tmp && mv .gitattributes.tmp .gitattributes

# Prompt for commit message
echo "Enter commit message (press Enter to use default message):"
read commit_message

# Use default message if no input provided
if [ -z "$commit_message" ]; then
  commit_message="Update .gitattributes and .gitignore for large files"
fi

# Commit changes to .gitattributes and .gitignore
git add .gitattributes .gitignore
git commit -m "$commit_message"

echo "Changes committed with message: $commit_message"
