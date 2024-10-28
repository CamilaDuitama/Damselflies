#!/bin/bash

# Set the size threshold (in bytes) for LFS tracking (e.g., 50MB)
SIZE_THRESHOLD=$((50 * 1024 * 1024))

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

# Commit changes to .gitattributes and .gitignore
git add .gitattributes .gitignore
git commit -m "Update .gitattributes and .gitignore for large files"
