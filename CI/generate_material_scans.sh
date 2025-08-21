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

# Function to create empty placeholder histograms for failed scans
create_empty_histograms() {
    local output_dir="$1"
    local compact_name="$2"
    local error_message="$3"
    
    echo "Creating empty placeholder histograms for failed scan: $compact_name"
    
    # Create ROOT macro to generate empty histograms
    cat > /tmp/create_empty_histograms.C << EOF
{
    // Create empty histograms for x0, lambda, and depth
    TString hist_types[] = {"x0", "lambda", "depth"};
    TString titles[] = {"X0 Material Budget [%]", "Lambda Interaction Length", "Material Depth"};
    
    for (int i = 0; i < 3; i++) {
        TString filename = "$output_dir/" + hist_types[i] + ".root";
        TFile *outfile = TFile::Open(filename, "RECREATE");
        
        if (!outfile || outfile->IsZombie()) {
            cout << "ERROR: Cannot create empty histogram file " << filename << endl;
            continue;
        }
        
        // Create empty histogram
        TH1F *empty_hist = new TH1F("empty_hist", "FAILED SCAN: $error_message", 100, 0, 180);
        empty_hist->SetTitle(titles[i] + " - SCAN FAILED: $compact_name");
        
        // Create canvas and stack for compatibility with existing pipeline
        TCanvas *canvas = new TCanvas("canvas", "Material Budget", 800, 600);
        THStack *stack = new THStack("stack", titles[i]);
        
        // Add the empty histogram to the stack
        empty_hist->SetFillColor(kRed);
        empty_hist->SetLineColor(kRed);
        stack->Add(empty_hist);
        
        // Draw and save
        canvas->cd();
        stack->Draw("hist");
        
        // Add text indicating failure
        TText *error_text = new TText(0.5, 0.5, "SCAN FAILED: $error_message");
        error_text->SetNDC();
        error_text->SetTextAlign(22);
        error_text->SetTextColor(kRed);
        error_text->SetTextSize(0.05);
        error_text->Draw();
        
        canvas->Write();
        outfile->Close();
        
        cout << "Created empty histogram: " << filename << endl;
    }
    
    cout << "Empty histograms created for failed scan: $compact_name" << endl;
}
EOF

    root -l -b -q /tmp/create_empty_histograms.C 2>&1
    rm -f /tmp/create_empty_histograms.C
}

# Function to process geometries and generate material scans directly into consolidated files
process_geometries() {
  local source_dir="$1"
  local consolidated_file="$2"
  local file_suffix="$3"
  
  echo "=== Processing geometries from $source_dir directly into $consolidated_file ==="
  
  # Initialize error tracking
  local failed_scans=()
  local total_processed=0
  local successful_scans=0
  
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
            total_processed=$((total_processed + 1))
            
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
            
            # Run material scan with error handling and timeout
            SCAN_SUCCESS=false
            
            echo "🔍 Starting material scan for ${compact_name}..."
            
            if [ "$QUIET_MODE" = true ]; then
              echo "Running material scan: $xml_file (k4run output suppressed)"
              
              # Use timeout and trap signals to catch segfaults
              timeout 300s bash -c "
                k4run utils/material_scan.py \
                  --GeoSvc.detector '$xml_file' \
                  --GeoDump.filename '${output_dir}/out_material_scan${file_suffix}.root' \
                  --angleDef theta \
                  --angleBinning $ANGLE_BINNING \
                  --angleMin $ANGLE_MIN \
                  --angleMax $ANGLE_MAX \
                  --nPhi $NPHI > /dev/null 2>&1
              " && SCAN_SUCCESS=true
              
            else
              echo "Running material scan: $xml_file"
              
              # Use timeout and trap signals to catch segfaults
              timeout 300s bash -c "
                k4run utils/material_scan.py \
                  --GeoSvc.detector '$xml_file' \
                  --GeoDump.filename '${output_dir}/out_material_scan${file_suffix}.root' \
                  --angleDef theta \
                  --angleBinning $ANGLE_BINNING \
                  --angleMin $ANGLE_MIN \
                  --angleMax $ANGLE_MAX \
                  --nPhi $NPHI
              " && SCAN_SUCCESS=true
              
            fi
            
            SCAN_EXIT_CODE=$?
            
            # Check if material scan succeeded
            if [ "$SCAN_SUCCESS" = true ] && [ -f "${output_dir}/out_material_scan${file_suffix}.root" ]; then
              echo "✅ Material scan completed successfully for ${compact_name}"
              
              # Generate plots normally
              echo "Generating material plots..."
              
              PLOT_SUCCESS=false
              timeout 120s python utils/material_plots.py \
                -f "${output_dir}/out_material_scan${file_suffix}.root" \
                -o "${output_dir}" \
                --angleDef theta \
                --angleBinning $ANGLE_BINNING \
                --angleMin $ANGLE_MIN \
                --angleMax $ANGLE_MAX && PLOT_SUCCESS=true
              
              if [ "$PLOT_SUCCESS" = true ] && [ -f "${output_dir}/x0.root" ]; then
                echo "✅ Plot generation completed successfully for ${compact_name}"
                successful_scans=$((successful_scans + 1))
              else
                echo "❌ Plot generation FAILED for ${compact_name}"
                failed_scans+=("${compact_name} (Plot generation failed)")
                create_empty_histograms "$output_dir" "$compact_name" "Plot generation failed"
              fi
              
            else
              # Material scan failed - determine why
              ERROR_REASON=""
              case $SCAN_EXIT_CODE in
                124) ERROR_REASON="Timeout (>5min)" ;;
                139) ERROR_REASON="Segmentation fault" ;;
                134) ERROR_REASON="SIGABRT signal" ;;
                137) ERROR_REASON="SIGKILL signal" ;;
                *) ERROR_REASON="Exit code $SCAN_EXIT_CODE" ;;
              esac
              
              echo "❌ MATERIAL SCAN FAILED for ${compact_name}: $ERROR_REASON"
              echo "🚨 ERROR: Material scan failure detected!" >&2
              
              failed_scans+=("${compact_name} ($ERROR_REASON)")
              
              # Create empty placeholder histograms
              create_empty_histograms "$output_dir" "$compact_name" "$ERROR_REASON"
            fi

            # Add histograms to consolidated file (empty or real)
            echo "Adding histograms to consolidated file for ${compact_name}"
            ./CI/combine_material_histograms.sh --add-to-consolidated \
              --input-dir "${output_dir}" \
              --output-file "$consolidated_file" \
              --histogram-prefix "${compact_name}" || {
              echo "❌ Failed to add histograms to consolidated file for ${compact_name}"
            }
            
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
  
  # Report final statistics
  echo ""
  echo "=== MATERIAL SCAN SUMMARY ==="
  echo "Total geometries processed: $total_processed"
  echo "Successful scans: $successful_scans"
  echo "Failed scans: ${#failed_scans[@]}"
  
  if [ ${#failed_scans[@]} -gt 0 ]; then
    echo ""
    echo "🚨 FAILED SCANS DETECTED:"
    for failed in "${failed_scans[@]}"; do
      echo "  ❌ $failed"
    done
    echo ""
    
    # Create error summary file for artifacts
    echo "# Material Scan Error Summary" > material_scan_errors.md
    echo "Date: $(date)" >> material_scan_errors.md
    echo "" >> material_scan_errors.md
    echo "## Summary" >> material_scan_errors.md
    echo "- Total processed: $total_processed" >> material_scan_errors.md
    echo "- Successful: $successful_scans" >> material_scan_errors.md
    echo "- Failed: ${#failed_scans[@]}" >> material_scan_errors.md
    echo "" >> material_scan_errors.md
    echo "## Failed Scans" >> material_scan_errors.md
    for failed in "${failed_scans[@]}"; do
      echo "- $failed" >> material_scan_errors.md
    done
    
    echo "📝 Error summary saved to: material_scan_errors.md"
  else
    echo "✅ All scans completed successfully!"
  fi
  echo "=========================="
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