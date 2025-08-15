#!/usr/bin/env python3

import argparse
import base64
import hashlib
import re
import sys


def hash_to_letters_md5(input_string):
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"  # uppercase only
    base = len(alphabet)

    # Step 1: Get MD5 hash as bytes
    hash_bytes = hashlib.md5(input_string.encode()).digest()

    # Step 2: Convert bytes to integer
    num = int.from_bytes(hash_bytes, 'big')

    # Step 3: Convert to base-52
    chars = []
    while num > 0:
        num, rem = divmod(num, base)
        chars.append(alphabet[rem])

    return ''.join(reversed(chars))[:8]

def process_taxonomy_string(taxonomy_string, prefix):
    """
    Process a taxonomy string and rename novel entries.
    
    Args:
        taxonomy_string (str): Original taxonomy string with levels separated by ';'
        prefix (str): Prefix to use for renamed entries
        
    Returns:
        str: Modified taxonomy string
    """
    levels = taxonomy_string.split(';')
    processed_levels = []
    
    for i, level in enumerate(levels):
        level = level.strip()
        
        # Check if this level matches the novel pattern
        if re.match(r'^[a-z]__novel_', level):
            # Extract the level character (first letter)
            level_char = level[0]
            
            # Build the taxonomy string up to this level (inclusive)
            taxonomy_up_to_here = ';'.join(processed_levels + [level])
            
            # Generate the base64 hash
            hash_suffix = hash_to_letters_md5(taxonomy_up_to_here)
            
            # Create the new level name
            new_level = f"{level_char}__{prefix}_{hash_suffix}"
            processed_levels.append(new_level)
        else:
            processed_levels.append(level)
    
    return ';'.join(processed_levels)


def main():
    parser = argparse.ArgumentParser(
        description='Rename novel taxonomy entries in vcontact taxonomy files'
    )
    parser.add_argument(
        '--vcontact-taxonomy', 
        required=True,
        help='Input TSV file with sequence names and vcontact taxonomy'
    )
    parser.add_argument(
        '--prefix', 
        required=True,
        help='Prefix string to use for renamed novel entries'
    )
    parser.add_argument(
        '--output-taxonomy',
        required=True,
        help='Output TSV file with modified taxonomy strings'
    )
    
    args = parser.parse_args()
    
    try:
        with open(args.vcontact_taxonomy, 'r') as infile, \
             open(args.output_taxonomy, 'w') as outfile:
            
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line:
                    continue
                    
                # Split into two columns
                parts = line.split('\t')
                if len(parts) != 2:
                    print(f"Warning: Line {line_num} does not have exactly 2 columns, skipping: {line}", 
                          file=sys.stderr)
                    continue
                
                sequence_name, taxonomy_string = parts
                
                # Process the taxonomy string
                modified_taxonomy = process_taxonomy_string(taxonomy_string, args.prefix)
                
                # Write to output
                outfile.write(f"{sequence_name}\t{modified_taxonomy}\n")
                
    except FileNotFoundError as e:
        print(f"Error: Could not find input file: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
