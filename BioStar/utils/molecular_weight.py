def get_weight_and_unit(weight: float) -> dict:
    if weight >= 1_000_000:
        transformed_weight = weight / 1_000_000.00
        unit = 'MDa'
    elif weight >= 1_000:
        transformed_weight = weight / 1_000.00
        unit = 'kDa'
    else:
        transformed_weight = weight
        unit = 'Da'

    return {
        'value': f"{transformed_weight:.2f}",
        'unit': unit
    }
