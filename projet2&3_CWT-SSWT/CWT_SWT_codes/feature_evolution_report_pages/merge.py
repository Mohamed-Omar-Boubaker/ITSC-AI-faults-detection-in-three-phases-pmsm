import os
from pypdf import PdfWriter

def merge_pdfs_in_folder(folder_path, output_filepath):
    merger = PdfWriter()
    
    # Retrieve and sort all PDF files in the directory
    pdf_files = sorted([f for f in os.listdir(folder_path) if f.lower().endswith('.pdf')])
    
    if not pdf_files:
        print("No PDF files found in the folder.")
        return

    # Append each PDF to the merger
    for file_name in pdf_files:
        file_path = os.path.join(folder_path, file_name)
        merger.append(file_path)
        print(f"Added: {file_name}")

    # Save the combined PDF
    merger.write(output_filepath)
    merger.close()
    print(f"\nSuccessfully created merged file: {output_filepath}")

if __name__ == "__main__":
    # Update these paths to match your folder and desired output location
    target_folder = "C:\Users\moham\OneDrive\Desktop\stage\codes\SWT_analyser_feats_extractor\feature_evolution_report_pages"
    output_pdf = "merged_output.pdf"
    
    merge_pdfs_in_folder(target_folder, output_pdf)