import numpy as np

def normalize_coords(*coords):
    """Normalized coordinates so that they are an array of (N, 2)"""

    if len(coords) == 1:
        if hasattr(coords[0], 'get_coords'):        # Catalog etc.
            return '', np.stack(coords[0].get_coords(), axis=-1)
        elif isinstance(coords[0], np.ndarray):
            return 'a', coords[0]
        elif isinstance(coords[0], tuple):
            return 't', stack_coords(coords[0][0], coords[0][1])
        elif isinstance(coords[0], list):
            return 'l', stack_coords(coords[0][0], coords[0][1])
        elif isinstance(coords[0], dict) and 'RA' in coords[0] and 'Dec' in coords[0]:
            return 'd', stack_coords(coords[0]['RA'], coords[0]['Dec'])
        elif isinstance(coords[0], dict) and 'x' in coords[0] and 'y' in coords[0]:
            return 'd', stack_coords(coords[0]['x'], coords[0]['y'])
        else:
            raise NotImplementedError('Coordinates cannot be interpreted as (ra, dec).')
    elif len(coords) == 2:
        return '2', np.stack(coords, axis=-1)
    else:
        raise ValueError('Coords must have 2 dimensions.')

def denormalize_coords(ctype, *coords):
    """Denormalizes coordinates from an array of (N, 2) into the format requested in `ctype`"""

    _, coords = normalize_coords(*coords)

    if ctype == '':
        return coords
    elif ctype == 'a':
        return coords
    elif ctype == 't':
        return (coords[..., 0], coords[..., 1])
    elif ctype == 'l':
        return [coords[..., 0], coords[..., 1]]
    elif ctype == 'd':
        return coords
    elif ctype == '2':
        return coords[..., 0], coords[..., 1]

    raise NotImplementedError()

def stack_coords(x, y):
    return np.stack([x, y], axis=-1)

def split_coords(*coords):
    ctype, coords = normalize_coords(*coords)
    return coords[..., 0], coords[..., 1]