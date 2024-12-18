import pandas as pd
import os
import sys
import argparse

def parse_arguments():
    """
    Parse arguments, including path, sample, and levels
    """
    parser = argparse.ArgumentParser(description = 'Extracts data at the specified classification level')
    parser.add_argument('-p', '--path', type = str, required = True,
                       help = 'The input & output path')
    parser.add_argument('-s', '--sample', type = str, required = True,
                       help = 'The name of the sample')
    parser.add_argument('-l', '--levels', nargs = '+', type = str,
                       help = 'The taxonomy levels to extract, could enter multiple, such as -l d p g')
    return parser.parse_args()

def validate_taxonomy_levels(input_levels):
    """
    Validate the inputted taxonomy levels
    """
    valid_levels = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    invalid_levels = [level for level in input_levels if level not in valid_levels]
    
    if invalid_levels:
        raise ValueError(f"Input taxonomy '{level}' is invalid. Valid taxonomies include {', '.join(valid_levels)}")
    
    if len(input_levels) > len(valid_levels):
        raise ValueError(f"Input taxonomy number ({len(input_levels)}) exceeded the maximum number of tanonomies ({len(valid_levels)})")
    
    return True

def get_taxonomy_levels(df):
    """
    Get all the unique prefixes
    """
    prefixes = []
    
    # Traverse each category string
    for taxonomy in df['Taxonomy']:
        # Separate different levels with |
        levels = taxonomy.split('|')
        
        # Extract the prefix (first letter) for each level
        for level in levels:
            prefix = level.strip().split('_')[0]
            if prefix not in prefixes:
                prefixes.append(prefix)
    
    return prefixes

def parse_taxonomy(taxonomy_string):
    """
    Parses the classification string
    Returning a dictionary with the key of the hierarchical type and the value of the corresponding classification name
    e.g. 'd__Bacteria|p__Proteobacteria' returns {'d': 'Bacteria', 'p': 'Proteobacteria'}
    """
    levels = taxonomy_string.split('|')
    result = {}
    for level in levels:
        level = level.strip()
        if '__' in level:
            prefix, name = level.split('__', 1)
            result[prefix] = name
    return result

def extract_taxonomy_level(df, target_level, count_col = 'Count', taxonomy_col = 'Taxonomy'):
    """
    Extract classification information and count data for a specified level
    
    Params:
    df: DataFrame, a data frame that contains classified data
    target_level: str, the target-level prefix (such as'd' , ' P' , ' C' , ' G' , 's' , etc.)
    count_col: str, the column name of the count data
    taxonomy_col: str, the column name of the taxonomy data
    
    Returns:
    DataFrame, containing the classification name and corresponding count data for the specified level
    """
    # Create the results list
    results = []
    
    # Traverse through each row of the data box
    for _, row in df.iterrows():
        taxonomy_dict = parse_taxonomy(row[taxonomy_col])
        
        # If a target level exists, it is added to the result
        if target_level in taxonomy_dict:
            results.append({
                'Level': target_level,
                'Name': taxonomy_dict[target_level],
                'Count': row[count_col]
            })
    
    # Convert to DataFrame and sort by count in descending order
    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df = result_df.sort_values('Count', ascending=False)
    
    return result_df

def save_level_matrix(df, output_file, target_level):
    """
    Saves the specified level of data to a file
    """
    level_df = extract_taxonomy_level(df, target_level)
    if not level_df.empty:
        with open(output_file, 'w') as f:
            for _, row in level_df.iterrows():
                f.write(f"{target_level}__{row['Name']}\t{row['Count']}\n")

def validate_taxonomy_level(level):
    valid_levels = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    if level not in valid_levels:
        raise ValueError(f"Input taxonomy '{level}' is invalid. Valid taxonomies include {', '.join(valid_levels)}")
    return True

if __name__ == "__main__":
    args = parse_arguments()

    path = args.path
    sample = args.sample

    dir = os.listdir( path )
    os.chdir( path )
    mpa_path = path + '/' + sample + '.mpa'
    mpa = pd.read_table(mpa_path, header = None, names = ['Taxonomy', 'Count'])

    # Extract bacteria taxonomy
    mpa_sub = mpa[ mpa.Taxonomy.str.contains("Bacteria") ]
    
    try:
        # Validate inputted taxonomy level(s)
        validate_taxonomy_levels(args.levels)
        
        # Generate output matrix with given taxonomy level(s)
        for level in args.levels:
            save_level_matrix(mpa_sub, f"counts_{level}.txt", level)
            print(f"Finish extracting {level} taxonomy read counts")
            
    except ValueError as e:
        print(f"Error: {e}")