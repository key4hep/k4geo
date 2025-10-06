#!/bin/bash

set -e  # Exit on any error

# Default values
QUIET_MODE=false
FAST_PARAMS=false
GEOMETRY_CONFIG=""

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -q, --quiet          Suppress verbose output from k4run material scanning"
    echo "  -f, --fast           Fast parameters: reduced angular resolution for quicker scanning"
    echo "  -c, --config FILE    Geometry configuration file (default: CI/config/geometry_list.yml)"
    echo "  -h, --help           Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                               # Use default geometry list"
    echo "  $0 --config CI/config/test.yml  # Use specific geometry list for testing"
    echo "  $0 --fast --quiet               # Fast scanning with default geometries"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -q|--quiet)
            QUIET_MODE=true
            shift
            ;;
        -f|--fast)
            FAST_PARAMS=true
            shift
            ;;
        -c|--config)
            GEOMETRY_CONFIG="$2"
            shift 2
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

# Set default config if not specified
if [ -z "$GEOMETRY_CONFIG" ]; then
    GEOMETRY_CONFIG="CI/config/geometry_list.yml"
fi

# Function to read geometry list from YAML config
read_geometry_list() {
    local config_file="$1"
    
    if [ ! -f "$config_file" ]; then
        echo "ERROR: Geometry configuration file not found: $config_file" >&2
        exit 1
    fi
    
    echo "Reading geometry list from: $config_file" >&2
    
    # Extract geometry paths from YAML (simple parsing)
    # This extracts lines that start with "  - " and removes the prefix
    local geometries=($(grep "^  - " "$config_file" | grep -v "^  - #" | sed 's/^  - //'))
    
    if [ ${#geometries[@]} -eq 0 ]; then
        echo "ERROR: No geometries found in config file: $config_file" >&2
        exit 1
    fi
    
    echo "Found ${#geometries[@]} geometries in config:" >&2
    for geom in "${geometries[@]}"; do
        echo "  - $geom" >&2
    done
    
    echo "${geometries[@]}"
}

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

# Modified function to process only configured geometries
process_geometries() {
  local source_dir="$1"
  local consolidated_file="$2"
  local file_suffix="$3"
  
  echo "=== Processing geometries from $source_dir directly into $consolidated_file ==="
  echo "Using configuration: $GEOMETRY_CONFIG"
  
  # Initialize error tracking
  local failed_scans=()
  local total_processed=0
  local successful_scans=0
  
  if [ "$FAST_PARAMS" = true ]; then
    echo "FAST MODE: Using reduced angular resolution"
  fi
  
  # Read geometry list from config
  local geometry_list=($(read_geometry_list "$GEOMETRY_CONFIG"))
  
  # Initialize the consolidated ROOT file
  ./CI/combine_material_histograms.sh --initialize "$consolidated_file"
  
  # Process each geometry from the config list
  for geometry_path in "${geometry_list[@]}"; do
    
    # Parse geometry path (format: GEOMETRY/COMPACT_NAME)
    local geometry_name=$(dirname "$geometry_path")
    local compact_name=$(basename "$geometry_path")
    
    # Construct full paths
    local geometry_dir="${source_dir}/FCCee/${geometry_name}"
    local compact_dir="${geometry_dir}/compact/${compact_name}"
    local xml_file="${compact_dir}/${compact_name}.xml"
    
    echo "Processing configured geometry: $geometry_path"
    
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
      echo "Warning: XML file not found for configured geometry: $xml_file"
      failed_scans+=("${compact_name} (XML file not found)")
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

# Function to detect the comparison target
detect_comparison_target() {
    local target_branch=""
    local target_repo=""
    
    echo "=== Detecting comparison target ==="
    echo "DEBUG: Environment variables:"
    echo "  GITHUB_BASE_REF='$GITHUB_BASE_REF'"
    echo "  GITHUB_REPOSITORY='$GITHUB_REPOSITORY'"
    echo "  GITHUB_HEAD_REF='$GITHUB_HEAD_REF'"
    echo "  PWD='$PWD'"
    
    # Check if we're in a git repository
    if git rev-parse --verify HEAD >/dev/null 2>&1; then
        echo "DEBUG: Git repository detected"
        echo "  Current branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'unknown')"
        echo "  Remote origin: $(git remote get-url origin 2>/dev/null || echo 'unknown')"
        echo "  Default branch: $(git symbolic-ref refs/remotes/origin/HEAD 2>/dev/null | sed 's@^refs/remotes/origin/@@' || echo 'unknown')"
    else
        echo "DEBUG: Not in a git repository"
    fi
    
    # Method 1: GitHub Actions environment (most reliable)
    if [ -n "$GITHUB_BASE_REF" ] && [ -n "$GITHUB_REPOSITORY" ]; then
        target_branch="$GITHUB_BASE_REF"
        target_repo="https://github.com/${GITHUB_REPOSITORY}.git"
        echo "✅ Detected from GitHub Actions:"
        echo "   Branch: $target_branch"
        echo "   Repository: $target_repo"
    
    # Method 2: Git commands (enhanced fallback)
    elif git rev-parse --verify HEAD >/dev/null 2>&1; then
        echo "Attempting to determine target branch from git..."
        
        # Try to get remote information
        if git remote get-url origin >/dev/null 2>&1; then
            target_repo=$(git remote get-url origin)
            echo "Found remote repository: $target_repo"
            
            # Try to fetch remote branch info
            if git ls-remote --symref origin HEAD >/dev/null 2>&1; then
                target_branch=$(git ls-remote --symref origin HEAD | head -1 | sed 's/^ref: refs\/heads\///' | cut -f1)
                echo "Found default branch: $target_branch"
            else
                echo "Could not determine default branch, using 'main'"
                target_branch="main"
            fi
        else
            echo "No git remote found, using defaults"
            target_branch="main"
            target_repo="https://github.com/key4hep/k4geo.git"
        fi
        
        echo "✅ Detected from git:"
        echo "   Branch: $target_branch"
        echo "   Repository: $target_repo"
    
    # Method 3: Defaults
    else
        target_branch="main"
        target_repo="https://github.com/key4hep/k4geo.git"
        echo "⚠️  Using defaults:"
        echo "   Branch: $target_branch"
        echo "   Repository: $target_repo"
    fi
    
    # Export for use in main script
    TARGET_BRANCH="$target_branch"
    TARGET_REPO="$target_repo"
}

# Call the detection function
detect_comparison_target

# Main execution
echo "=== Starting material histogram generation ==="

# Clone the main branch for reference data generation
echo "=== Cloning comparison target ==="
echo "Repository: $TARGET_REPO"
echo "Branch: $TARGET_BRANCH"

if [ "$QUIET_MODE" = true ]; then
  git clone --branch "$TARGET_BRANCH" --depth 1 "$TARGET_REPO" k4geo_main_ref > /dev/null 2>&1 || {
    echo "❌ Failed to clone $TARGET_BRANCH from $TARGET_REPO"
    echo "   Falling back to main branch from upstream..."
    git clone --branch main --depth 1 "https://github.com/key4hep/k4geo.git" k4geo_main_ref > /dev/null 2>&1 || {
      echo "❌ FATAL: Could not clone any reference branch"
      exit 1
    }
  }
else
  git clone --branch "$TARGET_BRANCH" --depth 1 "$TARGET_REPO" k4geo_main_ref || {
    echo "❌ Failed to clone $TARGET_BRANCH from $TARGET_REPO"
    echo "   Falling back to main branch from upstream..."
    git clone --branch main --depth 1 "https://github.com/key4hep/k4geo.git" k4geo_main_ref || {
      echo "❌ FATAL: Could not clone any reference branch"
      exit 1
    }
  }
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