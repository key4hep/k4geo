#!/bin/bash
# filepath: /home/fshaw/Documents/Code/k4geo/CI/gen_material_scan_hists.sh

# Default values
INPUT_DIR="."
OUTPUT_DIR="."

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -i, --input-dir DIR   Input directory containing ROOT files (default: current directory)"
    echo "  -o, --output-dir DIR  Output directory for summed histogram files (default: current directory)"
    echo "  -h, --help           Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -i /path/to/root/files -o /path/to/output"
    echo "  $0 --input-dir ./data --output-dir ./results"
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