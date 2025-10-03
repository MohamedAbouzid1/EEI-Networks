import numpy as np

def find_expression_with_mapping(coordinate, expression_df, coord_to_exon):
    """
    Find expression data using pre-loaded mappings.
    """
    # Guard against non-string or missing coordinates
    if coordinate is None or (isinstance(coordinate, float) and np.isnan(coordinate)):
        return None
    coordinate = str(coordinate).strip()
    if coordinate == '' or coordinate.lower() == 'nan':
        return None

    # Direct ENSMUSE ID lookup
    if coordinate.startswith('ENSMUSE'):
        if coordinate in expression_df.index:
            return expression_df.loc[coordinate]
        return None
    
    # Use mapping dictionary
    if coord_to_exon and coordinate in coord_to_exon:
        exon_id = coord_to_exon[coordinate]
        if exon_id in expression_df.index:
            return expression_df.loc[exon_id]
    
    # Try alternative coordinate formats
    if '-' in coordinate and coord_to_exon:
        # Convert chr1:12345-67890 to chr1:12345:67890
        parts = coordinate.replace('-', ':').split(':')
        if len(parts) == 3:
            # Try both strands
            for strand in ['1', '-1']:
                alt_coord = f"{parts[0]}:{parts[1]}:{parts[2]}:{strand}"
                if alt_coord in coord_to_exon:
                    exon_id = coord_to_exon[alt_coord]
                    if exon_id in expression_df.index:
                        return expression_df.loc[exon_id]
    
    return None