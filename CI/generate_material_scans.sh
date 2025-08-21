#!/bin/bash

set -e  # Exit on any error

# Default values
QUIET_MODE=false
TEST_MODE=false
FAST_PARAMS=false

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -q, --quiet          Suppress verbose output from k4run material scanning"
    echo "  -t, --test           Test mode: process only first 3 geometries and first compact file per geometry"
    echo "  -f, --fast           Fast parameters: reduced angular resolution for quicker scanning"
    echo "  -h, --help           Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                   # Full processing with verbose output"
    echo "  $0 --quiet           # Full processing with minimal output"
    echo "  $0 --test --quiet    # Test first 3 geometries with minimal output"
    echo "  $0 --fast --quiet    # Fast parameters with minimal output"
    echo "  $0 --test --fast     # Test mode with fast parameters (fastest)"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -q|--quiet)
            QUIET_MODE=true
            shift
            ;;
        -t|--test)
            TEST_MODE=true
            shift
            ;;
        -f|--fast)
            FAST_PARAMS=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Function to process geometries and generate material scans directly into consolidated files
process_geometries() {
  local source_dir="$1"
  local consolidated_file="$2"
  local file_suffix="$3"
  
  echo "=== Processing geometries from $source_dir directly into $consolidated_file ==="
  
  if [ "$TEST_MODE" = true ]; then
    echo "TEST MODE: Processing only first 3 geometries and first compact file per geometry"
  fi
  
  if [ "$FAST_PARAMS" = true ]; then
    echo "FAST MODE: Using reduced angular resolution"
  fi
  
  # Initialize the consolidated ROOT file
  ./CI/combine_material_histograms.sh --initialize "$consolidated_file"
  
  local geometry_count=0
  
  for geometry_dir in "${source_dir}"/FCCee/*/; do
    if [ -d "$geometry_dir" ]; then
      geometry_name=$(basename "$geometry_dir")
      echo "Processing geometry: $geometry_name"
      
      local compact_count=0
      
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

            output_dir="${source_dir}_temp/${geometry_name}/${compact_name}"
            mkdir -p "$output_dir"
            
            # Set parameters based on mode
            if [ "$FAST_PARAMS" = true ]; then
              # Fast test parameters - much reduced resolution
              ANGLE_BINNING=1
              ANGLE_MIN=0
              ANGLE_MAX=180
              NPHI=10
              echo "Using fast parameters: θ=${ANGLE_MIN}-${ANGLE_MAX}°, binning=${ANGLE_BINNING}°, φ=${NPHI}"
            else
              # Full parameters
              ANGLE_BINNING=1
              ANGLE_MIN=0
              ANGLE_MAX=180
              NPHI=100
            fi
            
            # Run material scan with optional output suppression
            if [ "$QUIET_MODE" = true ]; then
              echo "Running material scan: $xml_file (k4run output suppressed)"
              k4run utils/material_scan.py \
                --GeoSvc.detector "$xml_file" \
                --GeoDump.filename "${output_dir}/out_material_scan${file_suffix}.root" \
                --angleDef theta \
                --angleBinning $ANGLE_BINNING \
                --angleMin $ANGLE_MIN \
                --angleMax $ANGLE_MAX \
                --nPhi $NPHI > /dev/null 2>&1
            else
              echo "Running material scan: $xml_file"
              k4run utils/material_scan.py \
                --GeoSvc.detector "$xml_file" \
                --GeoDump.filename "${output_dir}/out_material_scan${file_suffix}.root" \
                --angleDef theta \
                --angleBinning $ANGLE_BINNING \
                --angleMin $ANGLE_MIN \
                --angleMax $ANGLE_MAX \
                --nPhi $NPHI
            fi
            
            # Generate plots directly into temporary directory
            echo "Generating material plots..."
            python utils/material_plots.py \
              -f "${output_dir}/out_material_scan${file_suffix}.root" \
              -o "${output_dir}" \
              --angleDef theta \
              --angleBinning $ANGLE_BINNING \
              --angleMin $ANGLE_MIN \
              --angleMax $ANGLE_MAX

            # Add histograms directly to consolidated file
            echo "Adding histograms to consolidated file for ${compact_name}"
            ./CI/combine_material_histograms.sh --add-to-consolidated \
              --input-dir "${output_dir}" \
              --output-file "$consolidated_file" \
              --histogram-prefix "${compact_name}"
            
            # Clean up temporary files
            rm -rf "${output_dir}"
            
            # In test mode, limit processing to first compact per geometry
            if [ "$TEST_MODE" = true ]; then
              compact_count=$((compact_count + 1))
              if [ $compact_count -ge 1 ]; then
                echo "Test mode: stopping after first compact directory for geometry $geometry_name"
                break
              fi
            fi
            
          else
            echo "Warning: XML file not found: $xml_file"
          fi
        fi
      done
      
      # In test mode, limit to first 3 geometries
      if [ "$TEST_MODE" = true ]; then
        geometry_count=$((geometry_count + 1))
        echo "Test mode: processed $geometry_count of 3 geometries"
        if [ $geometry_count -ge 3 ]; then
          echo "Test mode: stopping after first 3 geometries"
          break
        fi
      fi
    fi
  done
  
  # Finalize the consolidated ROOT file
  ./CI/combine_material_histograms.sh --finalize "$consolidated_file"
}

# Main execution
echo "=== Starting material histogram generation ==="

# Clone the main branch for reference data generation
echo "=== Cloning main branch for reference data ==="
if [ "$QUIET_MODE" = true ]; then
  git clone --branch implement-histcmp --depth 1 https://github.com/fredrikshaw/k4geo.git k4geo_main_ref > /dev/null 2>&1
else
  git clone --branch implement-histcmp --depth 1 https://github.com/fredrikshaw/k4geo.git k4geo_main_ref # THIS NEEDS TO BE CHANGED TO THE REAL MAIN BRANCH ON K4GEO BEFORE MERGE!!!!!!!
fi

# Generate reference histograms directly into consolidated file
cd k4geo_main_ref
process_geometries "." "../detector_geometries_ref.root" "_ref"
cd ..

# Generate current branch histograms directly into consolidated file
process_geometries "." "detector_geometries_monitored.root" ""

# Clean up temporary directories
rm -rf reference_results_temp mat_scan_results_temp k4geo_main_ref

echo "=== Material histogram generation completed ==="
echo "=== Consolidated files created: ==="
echo "  - detector_geometries_ref.root"
echo "  - detector_geometries_monitored.root"