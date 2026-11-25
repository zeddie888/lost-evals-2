import numpy as np
from scipy.spatial.transform import Rotation

def radec_roll_to_quat(ra_deg: float, de_deg: float, roll_deg: float) -> Rotation:
    """Build the inertial-to-camera-frame rotation for a given (ra, de, roll),
    matching LOST's internal convention (validated against real LOST output):
    boresight is axis 0 (not axis 2), and roll rotates (u_ref, v_ref) by the
    matrix below.
    """
    ra, de, roll = np.radians([ra_deg, de_deg, roll_deg])

    boresight = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])

    north = np.array([0.0, 0.0, 1.0])
    v_ref = north - np.dot(north, boresight) * boresight
    v_ref /= np.linalg.norm(v_ref)
    u_ref = np.cross(v_ref, boresight)

    row1 = np.cos(roll) * u_ref - np.sin(roll) * v_ref
    row2 = np.sin(roll) * u_ref + np.cos(roll) * v_ref

    R = np.array([boresight, row1, row2])
    return Rotation.from_matrix(R)

def cross_and_roll_error_arcsec(true_radec_roll, est_quat_ijkw) -> tuple[float, float]:
    """Decompose the rotation error between true attitude (ra, de, roll) and
    LOST's estimated attitude (quaternion, i/j/k/real order) into
    cross-boresight and around-boresight error, in arcseconds.
    """
    r_true = radec_roll_to_quat(*true_radec_roll)
    r_est = Rotation.from_quat(est_quat_ijkw)

    r_err_body = r_est * r_true.inv()
    rotvec_rad = r_err_body.as_rotvec()

    roll_arcsec = rotvec_rad[0] * 206265.0    # boresight
    pitch_arcsec = rotvec_rad[1] * 206265.0
    yaw_arcsec = rotvec_rad[2] * 206265.0

    return pitch_arcsec, yaw_arcsec, roll_arcsec