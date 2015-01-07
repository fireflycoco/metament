import numpy as np

__all__ = ['wswd2uv', 'uv2wswd']

# Using np.multiply and np.add and np.square to keep scaler.
def wswd2uv(ws, wd):
    """Convert wind speed and wind direction to u, v.
    ws: wind speed
    wd: wind direction (in degrees, north wind is 0, east wind is 90, etc)
    Returns: u, v
    """
    wd = np.deg2rad(wd)
    u = -np.multiply(ws, np.sin(wd))
    v = -np.multiply(ws, np.cos(wd))
    return u, v

def uv2wswd(u, v):
    """Convert u, v wind to wind speed and wind direction.
    u , v : u, v wind.
    Returns: ws, wd
    ws: wind speed
    wd: wind direction (in degrees, north wind is 0, east wind is 90, etc)
    """
    ws = np.hypot(u, v)
    wd = np.fmod(np.rad2deg(np.arctan2(u, v)) + 180.0 , 360.0)
    return ws, wd

