#!/bin/bash
# filepath: /home/fshaw/Documents/Code/k4geo/CI/generate_material_scans.sh

set -e  # Exit on any error

# Default values
QUIET_MODE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --quiet|-q)
            QUIET_MODE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Function to process geometries and generate material scans
process_geometries() {
  local source_dir="$1"
  local output_base="$2"
  local file_suffix="$3"
  
  echo "=== Processing geometries from $source_dir ==="
  
  for geometry_dir in "${source_dir}"/FCCee/*/; do
    if [ -d "$geometry_dir" ]; then
      geometry_name=$(basename "$geometry_dir")
      echo "Processing geometry: $geometry_name"
      
      # Iterate through all compact directories for this geometry
      for compact_dir in "${geometry_dir}compact"/*/; do
        if [ -d "$compact_dir" ]; then
          compact_name=$(basename "$compact_dir")
          xml_file="${compact_dir}${compact_name}.xml"
          
          # Check if the XML file exists
          if [ -f "$xml_file" ]; then
            if [ "$QUIET_MODE" = true ]; then
              echo "Processing: $xml_file (output suppressed)"
            else
              echo "Processing: $xml_file"
            fi

            # ADD THIS LINE - Define output directory
            output_dir="${output_base}/${geometry_name}/${compact_name}"
            mkdir -p "$output_dir"
            
            # Run material scan with optional output suppression
            if [ "$QUIET_MODE" = true ]; then
              echo "Running material scan: $xml_file (k4run output suppressed)"
              k4run utils/material_scan.py \
                --GeoSvc.detector "$xml_file" \
                --GeoDump.filename "${output_dir}/out_material_scan${file_suffix}.root" \
                --angleDef theta \
                --angleBinning 1 \
                --angleMin 0 \
                --angleMax 180 \
                --nPhi 100 > /dev/null 2>&1
            else
              echo "Running material scan: $xml_file"
              k4run utils/material_scan.py \
                --GeoSvc.detector "$xml_file" \
                --GeoDump.filename "${output_dir}/out_material_scan${file_suffix}.root" \
                --angleDef theta \
                --angleBinning 1 \
                --angleMin 0 \
                --angleMax 180 \
                --nPhi 100
            fi
            
            # Generate plots (always show output)
            echo "Generating material plots..."
            python utils/material_plots.py \
              -f "${output_dir}/out_material_scan${file_suffix}.root" \
              -o "${output_dir}" \
              --angleDef theta \
              --angleBinning 1 \
              --angleMin 0 \
              --angleMax 180

            # Rename plot files if this is reference data
            if [ "$file_suffix" = "_ref" ]; then
              cd "$output_dir"
              for file in x0.root lambda.root depth.root; do
                if [ -f "$file" ]; then
                  mv "$file" "${file%.root}${file_suffix}.root"
                fi
              done
              cd - > /dev/null
            fi

            # Combine histogram channels
            ./CI/combine_material_histograms.sh --input-dir "${output_dir}" --output-dir "${output_dir}"
            
            # Rename summed files if this is reference data
            if [ "$file_suffix" = "_ref" ]; then
              cd "$output_dir"
              for file in *_summed_hist.root; do
                if [ -f "$file" ]; then
                  mv "$file" "${file%_summed_hist.root}${file_suffix}_summed_hist.root"
                fi
              done
              cd - > /dev/null
            fi
            
          else
            echo "Warning: XML file not found: $xml_file"
          fi
        fi
      done
    fi
  done
}

# Main execution
echo "=== Starting material histogram generation ==="

# Clone the main branch for reference data generation
echo "=== Cloning main branch for reference data ==="
git clone --branch implement-histcmp --depth 1 https://github.com/fredrikshaw/k4geo.git k4geo_main_ref # THIS NEEDS TO BE CHANGED TO THE REAL MAIN BRANCH ON K4GEO BEFORE MERGE!!!!!!!

# Generate reference histograms from main branch
cd k4geo_main_ref
process_geometries "." "../reference_results" "_ref"
cd ..

# Generate current branch histograms
process_geometries "." "mat_scan_results" ""

# Copy reference files to current branch directories for comparison
echo "=== Copying reference files for comparison ==="
for geometry_dir in mat_scan_results/*/; do
  if [ -d "$geometry_dir" ]; then
    geometry_name=$(basename "$geometry_dir")
    
    for compact_dir in "${geometry_dir}"*/; do
      if [ -d "$compact_dir" ]; then
        compact_name=$(basename "$compact_dir")
        ref_source_dir="reference_results/${geometry_name}/${compact_name}"
        
        if [ -d "$ref_source_dir" ]; then
          echo "Copying reference files for ${geometry_name}/${compact_name}"
          cp "$ref_source_dir"/*_ref*.root "$compact_dir/" 2>/dev/null || echo "No reference files to copy"
        fi
      fi
    done
  fi
done

echo "=== Material histogram generation completed ==="