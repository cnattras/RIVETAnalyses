#!/bin/bash


# Define the filename within the script
filename="PHENIX_2008_I777211.plot"

# Check if the file exists
if [ ! -f "$filename" ]; then
  echo "Error: File '$filename' not found!"
  exit 1
fi

# Define a variable to indicate if the DrawOnly section is found
drawonly_section_found=false

# Read the file line by line
while IFS= read -r line; do
  # Check if the line contains the string "DrawOnly"
  if [[ "$line" == DoNotDraw:* ]]; then
    # Extract the text following "DrawOnly"
    drawonly_content="${line#DoNotDraw: }"
    # Output the content
    echo "DrawOnly content: $drawonly_content"
    # Set the flag to true
    drawonly_section_found=true
    # Exit the loop since we found the DrawOnly content
    break
  fi
done < "$filename"

# Check if the DrawOnly section was not found
if ! $drawonly_section_found; then
  rivet-mkhtml --pwd Rivet.yoda
fi

if $drawonly_section_found; then
  # Separate the content by whitespace
  read -ra drawonly_array <<< "$drawonly_content"
  # Define a variable to store the echoed content
  output=""
  # Populate the echoed content
  for element in "${drawonly_array[@]}"; do
    output+="-M $element "
    done
fi

rivet-mkhtml --pwd $output Rivet.yoda
