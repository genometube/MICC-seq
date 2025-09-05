#!/bin/bash

OUTPUT_FILE="project_documentation.md"
DIR_PATH="/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/test/wj_n3kJQ1_250829"

{
  echo "# Project Documentation"
  echo ""
  echo "## Text Content"
  echo ""
  
  # Process text files
  for txt_file in "$DIR_PATH"/*.txt; do
    echo "### $(basename "$txt_file")"
    echo "\`\`\`"
    cat "$txt_file"
    echo "\`\`\`"
    echo ""
  done
  
  # Process TSV files as markdown tables
  echo "## Data Metrics"
  echo ""
  for tsv_file in "$DIR_PATH"/*.tsv; do
    echo "### $(basename "$tsv_file")"
    echo ""
    
    # Get headers
    headers=$(head -n 1 "$tsv_file")
    echo "| $headers |" | sed 's/\t/ | /g'
    
    # Add separator line
    num_cols=$(echo "$headers" | awk -F'\t' '{print NF}')
    sep_line=$(printf "| %s |" $(seq -s " | " $num_cols | sed 's/[0-9]\+/---/g'))
    echo "$sep_line"
    
    # Add data rows
    tail -n +2 "$tsv_file" | while IFS=$'\t' read -r -a line; do
      echo "| ${line[*]} |" | sed 's/\t/ | /g'
    done
    echo ""
  done
  
  # Process PNG files
  echo "## Visualizations"
  echo ""
  count=1
  for png_file in "$DIR_PATH"/*.png; do
    echo "![Figure $count]($(basename "$png_file"))"
    echo "*Figure $count: $(basename "$png_file" | sed 's/\.[^.]*$//' | tr '_' ' ')*"
    echo ""
    ((count++))
  done
  
  echo "Generated on $(date)"
} > "$DIR_PATH/$OUTPUT_FILE"

echo "Markdown file created: $DIR_PATH/$OUTPUT_FILE"