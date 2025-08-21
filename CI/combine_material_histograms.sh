#!/bin/bash

# Default values
INPUT_DIR="."
OUTPUT_DIR="."
MODE="legacy"
OUTPUT_FILE=""
HISTOGRAM_PREFIX=""

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -i, --input-dir DIR          Input directory containing ROOT files (default: current directory)"
    echo "  -o, --output-dir DIR         Output directory for summed histogram files (default: current directory)"
    echo "  --initialize FILE            Initialize a new consolidated ROOT file"
    echo "  --add-to-consolidated        Add histograms to existing consolidated file"
    echo "  --output-file FILE           Target consolidated ROOT file"
    echo "  --histogram-prefix PREFIX    Prefix for histogram names"
    echo "  --finalize FILE              Finalize consolidated ROOT file"
    echo "  -h, --help                   Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 -i /path/to/root/files -o /path/to/output  # Legacy mode"
    echo "  $0 --initialize detector_geometries.root      # Initialize consolidated file"
    echo "  $0 --add-to-consolidated --input-dir ./data --output-file detector_geometries.root --histogram-prefix ALLEGRO_o1_v01"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --initialize)
            MODE="initialize"
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --add-to-consolidated)
            MODE="add"
            shift
            ;;
        --output-file)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --histogram-prefix)
            HISTOGRAM_PREFIX="$2"
            shift 2
            ;;
        --finalize)
            MODE="finalize"
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            show_usage
            exit 1
            ;;
    esac
done

# Function to initialize consolidated ROOT file
initialize_consolidated_file() {
    local output_file="$1"
    echo "Initializing consolidated ROOT file: $output_file"
    
    cat > /tmp/initialize_consolidated.C << EOF
{
    TFile *outfile = TFile::Open("$output_file", "RECREATE");
    if (!outfile || outfile->IsZombie()) {
        cout << "ERROR: Cannot create consolidated file $output_file" << endl;
        return;
    }
    
    cout << "Initialized consolidated file: $output_file" << endl;
    outfile->Close();
}
EOF

    root -l -b -q /tmp/initialize_consolidated.C 2>&1
    rm -f /tmp/initialize_consolidated.C
}

# Function to add histograms to consolidated file
add_to_consolidated_file() {
    local input_dir="$1"
    local output_file="$2" 
    local hist_prefix="$3"
    
    echo "Adding histograms from $input_dir to $output_file with prefix $hist_prefix"
    
    cat > /tmp/add_to_consolidated.C << EOF
{
    // Open consolidated file in UPDATE mode
    TFile *outfile = TFile::Open("$output_file", "UPDATE");
    if (!outfile || outfile->IsZombie()) {
        cout << "ERROR: Cannot open consolidated file $output_file" << endl;
        return;
    }
    
    cout << "Adding histograms with prefix: $hist_prefix" << endl;
    
    // Process each histogram type
    TString hist_types[] = {"x0", "lambda", "depth"};
    for (int i = 0; i < 3; i++) {
        TString input_file = "$input_dir/" + hist_types[i] + ".root";
        
        // Open input file
        TFile *infile = TFile::Open(input_file, "READ");
        if (!infile || infile->IsZombie()) {
            cout << "WARNING: Cannot open " << input_file << endl;
            continue;
        }
        
        // Find canvas and extract THStack
        TCanvas *canvas = nullptr;
        TIter next(infile->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next())) {
            TObject *obj = key->ReadObj();
            if (obj->IsA() == TCanvas::Class()) {
                canvas = (TCanvas*)obj;
                break;
            }
        }
        
        if (!canvas) {
            cout << "WARNING: No canvas found in " << input_file << endl;
            infile->Close();
            continue;
        }
        
        // Get THStack from canvas
        if (canvas->GetListOfPrimitives()->GetEntries() < 2) {
            cout << "WARNING: Canvas has insufficient primitives in " << input_file << endl;
            infile->Close();
            continue;
        }
        
        TObject *obj1 = canvas->GetListOfPrimitives()->At(1);
        THStack *stack = (THStack*)obj1;
        if (!stack || !stack->InheritsFrom("THStack")) {
            cout << "WARNING: No THStack found in " << input_file << endl;
            infile->Close();
            continue;
        }
        
        // Get histograms from stack and combine them
        TList *hists = stack->GetHists();
        if (!hists || hists->GetEntries() == 0) {
            cout << "WARNING: No histograms in stack from " << input_file << endl;
            infile->Close();
            continue;
        }
        
        // Create combined histogram
        TH1F *combined = (TH1F*)hists->At(0)->Clone();
        for (int j = 1; j < hists->GetEntries(); j++) {
            TH1F *hist = (TH1F*)hists->At(j);
            combined->Add(hist);
        }
        
        // Set name and title for consolidated file
        TString hist_name = "$hist_prefix-" + hist_types[i];
        combined->SetName(hist_name);
        combined->SetTitle("Material Budget: $hist_prefix " + hist_types[i]);
        
        // Write to consolidated file
        outfile->cd();
        combined->Write();
        cout << "Added histogram: " << hist_name << endl;
        
        infile->Close();
    }
    
    outfile->Close();
    cout << "Successfully added histograms for $hist_prefix" << endl;
}
EOF

    root -l -b -q /tmp/add_to_consolidated.C 2>&1
    rm -f /tmp/add_to_consolidated.C
}

# Function to finalize consolidated file
finalize_consolidated_file() {
    local output_file="$1"
    echo "Finalizing consolidated ROOT file: $output_file"
    
    cat > /tmp/finalize_consolidated.C << EOF
{
    TFile *file = TFile::Open("$output_file", "READ");
    if (!file || file->IsZombie()) {
        cout << "ERROR: Cannot open $output_file for verification" << endl;
        return;
    }
    
    cout << "=== Final consolidated file contents ===" << endl;
    file->ls();
    cout << "=========================================" << endl;
    
    file->Close();
    cout << "Consolidated file finalized: $output_file" << endl;
}
EOF

    root -l -b -q /tmp/finalize_consolidated.C 2>&1
    rm -f /tmp/finalize_consolidated.C
}

# Main execution based on mode
case $MODE in
    "initialize")
        initialize_consolidated_file "$OUTPUT_FILE"
        ;;
    "add")
        if [ -z "$OUTPUT_FILE" ] || [ -z "$HISTOGRAM_PREFIX" ]; then
            echo "Error: --output-file and --histogram-prefix required for add mode"
            exit 1
        fi
        add_to_consolidated_file "$INPUT_DIR" "$OUTPUT_FILE" "$HISTOGRAM_PREFIX"
        ;;
    "finalize")
        finalize_consolidated_file "$OUTPUT_FILE"
        ;;
    "legacy")
        # Keep existing functionality for backward compatibility
        
        # Validate directories exist
        if [ ! -d "$INPUT_DIR" ]; then
            echo "Error: Input directory $INPUT_DIR does not exist"
            exit 1
        fi

        # Create output directory if it doesn't exist
        if [ ! -d "$OUTPUT_DIR" ]; then
            echo "Creating output directory: $OUTPUT_DIR"
            mkdir -p "$OUTPUT_DIR"
        fi

        echo "Using input directory: $INPUT_DIR"
        ls -la "$INPUT_DIR" || echo "No files found in input directory"
        echo "Using output directory: $OUTPUT_DIR"

        # Function to process a ROOT file and create combined histogram
        process_root_file() {
            local input_file="$1"
            local output_file="$2"
            local hist_title="$3"
            
            # echo "Input file directory: $(dirname "$input_file")"
            # ls -la "$(dirname "$input_file")"

            if [ ! -f "$input_file" ]; then
                echo "Warning: $input_file not found, skipping..."
                return
            fi
            
            echo "Processing $input_file..."
            
            # Create temporary ROOT macro
            cat > /tmp/process_root.C << EOF
{
    // Open the input file
    cout << "Attempting to open: $input_file" << endl;
    TFile *infile = TFile::Open("$input_file", "READ");
    if (!infile || infile->IsZombie()) {
        cout << "ERROR: Cannot open $input_file" << endl;
        return;
    }

    cout << "File opened successfully. Listing contents:" << endl;
    infile->ls();

    // Get the canvas
    TCanvas *canvas = nullptr;
    cout << "Looking for canvas..." << endl;
    TIter next(infile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        cout << "Found object: " << obj->GetName() << " of class " << obj->ClassName() << endl;
        if (obj->IsA() == TCanvas::Class()) {
            canvas = (TCanvas*)obj;
            cout << "Found canvas: " << canvas->GetName() << endl;
            break;
        }
    }

    if (!canvas) {
        cout << "ERROR: No canvas found in $input_file" << endl;
        infile->Close();
        return;
    }

    // List canvas contents
    cout << "Canvas primitives:" << endl;
    canvas->GetListOfPrimitives()->ls();

    // Try to get THStack
    cout << "Attempting to get THStack at position 1..." << endl;
    if (canvas->GetListOfPrimitives()->GetEntries() < 2) {
        cout << "ERROR: Canvas has less than 2 primitives" << endl;
        infile->Close();
        return;
    }

    TObject *obj1 = canvas->GetListOfPrimitives()->At(1);
    cout << "Object at position 1: " << obj1->GetName() << " class: " << obj1->ClassName() << endl;

    THStack *stack = (THStack*)obj1;
    if (!stack || !stack->InheritsFrom("THStack")) {
        cout << "ERROR: Object at position 1 is not a THStack" << endl;
        infile->Close();
        return;
    }

    // Get histograms from stack
    TList *hists = stack->GetHists();
    if (!hists || hists->GetEntries() == 0) {
        cout << "ERROR: No histograms found in stack" << endl;
        infile->Close();
        return;
    }

    cout << "Found " << hists->GetEntries() << " histograms in stack" << endl;

    // Create combined histogram
    TH1F *combined = (TH1F*)hists->At(0)->Clone("combined_hist");
    combined->SetTitle("$hist_title");
    combined->SetName("combined_hist");

    // Sum all other histograms
    for (int i = 1; i < hists->GetEntries(); i++) {
        TH1F *hist = (TH1F*)hists->At(i);
        combined->Add(hist);
    }

    // Create output file
    cout << "Creating output file: $output_file" << endl;
    TFile *outfile = TFile::Open("$output_file", "RECREATE");
    if (!outfile || outfile->IsZombie()) {
        cout << "ERROR: Cannot create output file $output_file" << endl;
        infile->Close();
        return;
    }

    combined->Write();
    outfile->Close();

    cout << "SUCCESS: Created combined histogram with " << hists->GetEntries() << " histograms" << endl;
    cout << "Combined histogram entries: " << combined->GetEntries() << endl;
    cout << "Saved to: $output_file" << endl;

    infile->Close();
}
EOF

            # Execute the ROOT macro and capture output
            echo "Executing ROOT macro..."
            root -l -b -q /tmp/process_root.C 2>&1
            
            # Clean up
            rm -f /tmp/process_root.C
            
            echo "Finished processing $input_file"
            echo "----------------------------------------"
        }

        # Process each type of ROOT file
        echo "=== Generating combined material scan histograms ==="

        # Process x0.root
        process_root_file "$INPUT_DIR/x0.root" "$OUTPUT_DIR/x0_summed_hist.root" "Combined X0 Material Budget x/X [%]"

        # Process lambda.root  
        process_root_file "$INPUT_DIR/lambda.root" "$OUTPUT_DIR/lambda_summed_hist.root" "Combined Lambda Interaction Length"

        # Process depth.root
        process_root_file "$INPUT_DIR/depth.root" "$OUTPUT_DIR/depth_summed_hist.root" "Combined Depth Material Scan"

        echo "=== Finished processing material scan histograms ==="

        # Optional: List the created files
        echo "Created files:"
        ls -la "$OUTPUT_DIR"/*_summed_hist.root 2>/dev/null || echo "No summed histogram files were created"

        echo "=== Legacy mode processing completed ==="
        ;;
    *)
        echo "Error: Invalid mode"
        exit 1
        ;;
esac