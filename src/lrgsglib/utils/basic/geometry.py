from .common import *

def project_3d_to_2d(x, y, z, theta=0., phi=0.):
    """
    Projects a 3D point (x, y, z) onto a 2D plane using specified rotation angles.

    Parameters
    ----------
    x : float
        The x-coordinate of the 3D point.
    y : float
        The y-coordinate of the 3D point.
    z : float
        The z-coordinate of the 3D point.
    theta : float, optional
        Rotation angle around the y-axis, in radians. If not provided, uses self.theta.
    phi : float, optional
        Rotation angle around the x-axis, in radians. If not provided, uses self.phi.

    Returns
    -------
    tuple of float
        The (x, y) coordinates of the point projected onto the 2D plane.
    """
    # Rotation matrix around the y-axis (theta)
    R_theta = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

    # Rotation matrix around the x-axis (phi)
    R_phi = np.array([
        [1, 0, 0],
        [0, np.cos(phi), -np.sin(phi)],
        [0, np.sin(phi), np.cos(phi)]
    ])

    # Initial position vector
    position = np.array([x, y, z])

    # Apply rotations (order matters)
    position_rotated = R_phi @ R_theta @ position

    # Project onto 2D plane (ignore z after rotation)
    x2, y2 = position_rotated[0], position_rotated[1]

    return x2, y2