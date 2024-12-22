from numpy import ma, zeros, where


def check_and_remove_masked(table, logger=None):
    """
    Check for masked columns in the table and remove rows with masked values.
    Parameters:
    table (Table): An astropy Table object to be checked and cleaned.
    Returns:
    Table: A new Table object with rows containing masked values removed.
    Notes:
    - This function checks if any columns in the table are MaskedColumns.
    - If no masked columns are found, the original table is returned.
    - If masked columns are found, rows with masked values are removed.
    - The function prints the number of rows removed due to masked values.
    """

    # See if any columns in the table are MaskedColumns
    masked = False
    for col in table.columns:
        if ma.is_masked(table[col]):
            masked = True
            break
    if not masked:
        return table
    
    # Get indices of rows with masked values
    mask = zeros(len(table), dtype=bool)
    for col in table.columns:
        mask |= ma.getmaskarray(table[col])
    good_rows = where(mask == 0)[0]

    table = table[good_rows]
    for col in table.columns:
        table[col] = table[col].data
    

    if logger is not None:
        logger.info(f'Removed {len(mask) - len(good_rows)} rows with masked values')
    else:
        print(f'Removed {len(mask) - len(good_rows)} rows with masked values')

    return table
