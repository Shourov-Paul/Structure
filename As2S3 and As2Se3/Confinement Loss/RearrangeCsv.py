import os
import csv

def rearrange():
    input_path = r"d:\Structure\As2S3 and As2Se3\Confinement Loss\All_n_Confloss_R2.5(um).csv"
    output_path = r"d:\Structure\As2S3 and As2Se3\Confinement Loss\All_n_Confloss_R2.5(um)_rearranged.csv"
    
    with open(input_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    if not lines:
        print("Empty file")
        return
        
    header = lines[0].strip()
    data_lines = [line.strip() for line in lines[1:] if line.strip()]
    
    # Let's count columns. Each line is split by comma.
    # The header might have a trailing comma, so let's keep that structure.
    num_rows = len(data_lines)
    print(f"Number of data rows: {num_rows}")
    
    # We want to parse columns.
    # Since there are trailing commas, let's split each line.
    row_data = [line.split(',') for line in data_lines]
    num_cols = len(row_data[0])
    print(f"Number of columns (including empty trailing ones): {num_cols}")
    
    # For each column (except the last empty one), we want to extract non-empty values,
    # shift them to the bottom, and fill the top with "0".
    columns_data = [[] for _ in range(num_cols)]
    for row in row_data:
        for col_idx in range(num_cols):
            # Ensure we don't index out of bounds
            val = row[col_idx] if col_idx < len(row) else ""
            columns_data[col_idx].append(val)
            
    # Rearrange columns
    new_columns_data = [[] for _ in range(num_cols)]
    
    # For col 4 (which is always empty string), we keep it empty.
    # For others, we shift non-empty values to the bottom.
    for col_idx in range(num_cols):
        col_vals = columns_data[col_idx]
        
        # If it is the last column and all values are empty, keep it empty
        if col_idx == num_cols - 1 and all(v == "" for v in col_vals):
            new_columns_data[col_idx] = [""] * num_rows
            continue
            
        non_empty = [v for v in col_vals if v != ""]
        num_non_empty = len(non_empty)
        num_zeros = num_rows - num_non_empty
        
        new_col = ["0"] * num_zeros + non_empty
        new_columns_data[col_idx] = new_col
        
    # Now write the new CSV file
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        # Write header
        f.write(header + "\n")
        # Write data rows
        for row_idx in range(num_rows):
            row_vals = [new_columns_data[col_idx][row_idx] for col_idx in range(num_cols)]
            f.write(",".join(row_vals) + "\n")
            
    print("Rearrangement complete.")

if __name__ == "__main__":
    rearrange()

