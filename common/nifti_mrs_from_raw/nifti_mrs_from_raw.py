""" Conversion of Bruker rawdata.job0 to NIfTI-MRS

Copyright:
William Clarke and Clemence Ligneul
University of Oxford
2024
"""

from pathlib import Path
import re

import numpy as np

from nifti_mrs import create_nmrs


def _load_array(data_path, nfids=None, npoints=None, points_to_remove=76):
    """Read a raw.job0 file (and probably others)

    Either specify nfids (number of transients), or
    npoints (number of time/spectral points).

    :param data_path: Path to file to read
    :type data_path: str or pathlib.Path
    :param nfids: Number of FIDs expected
    :type nfids: int, optional
    :param nfids: Number of time points expected
    :type nfids: int, optional
    :param points_to_remove:
        Number of points to remove from the start of each FID, defaults to 76
    :type points_to_remove: int, optional
    :return: Array of complex FIDs
    :rtype: numpy.ndarray
    """
    raw_data = np.fromfile(data_path, dtype=np.int32).astype(float)
    raw_data_cmplx = raw_data[0::2] + 1j * raw_data[1::2]

    if nfids is None and npoints is None:
        raise ValueError('One of nfids or npoints must be specified.')
    elif nfids is None:
        nfids = int(raw_data_cmplx.shape[0] / npoints)
    elif npoints is None:
        npoints = int(raw_data_cmplx.shape[0] / nfids)
    raw_data_cmplx = raw_data_cmplx.reshape(nfids, npoints)

    return np.concatenate(
        (raw_data_cmplx[:, points_to_remove:],
         np.zeros((raw_data_cmplx.shape[0], points_to_remove))),
        axis=1)


def _read_header_file_info(file_path, keys_single, keys_array):
    """Read information from the method file

    :param file_path: path to the header file
    :type file_path: str or pathlib.Path
    :param keys_single: List of header keys that are a single value
    :type keys_single: list of str
    :param keys_array: List of header keys that have array values
    :type keys_array: list of str
    :return: Dict containing the information
    :rtype: dict
    """
    re_searches = [re.compile(fr'({x})\=(\d+)') for x in keys_single]
    re_searches2 = [re.compile(fr'({x})\=\((\s?\d+\s?)\)') for x in keys_array]

    with open(file_path) as fp:
        methodlines = fp.readlines()

    method_values = {}
    for line in methodlines:
        for re_ptrn in re_searches:
            match = re.search(re_ptrn, line)
            if match:
                method_values[match[1]] = int(match[2])

    # For array values that occur on the line after
    for idx, line in enumerate(methodlines):
        for re_ptrn in re_searches2:
            match = re.search(re_ptrn, line)
            if match:
                method_values[match[1]] = float(
                    methodlines[idx+1].split(' ')[0])

    return method_values


def _correct_offset(data, dwell, offset_hz):
    """Apply linear phase to correct a frequency offset

    :param data: data (nfids * npoints)
    :type data: np.ndarray
    :param dwell: dwell time in seconds
    :type dwell: float
    :param offset_hz: Frequency offset to correct, in hertz
    :type offset_hz: float
    :return: shifted data
    :rtype: np.ndarray
    """
    time_axis = np.arange(data.shape[1]) * dwell
    lin_phase = np.exp(1j * time_axis * offset_hz * 2 * np.pi)
    return data * lin_phase


def main(path, out, nfids=None, correct_work_offset=True):
    """Reads a rawdata.job0 and outputs a NIfTI-MRS formated file

    If the number of FIDs matches the stored data size, this function
    will separate PVM_EncNReceivers, PVM_NAverages, and PVM_NRepetitions
    into 'DIM_COIL', 'DIM_DYN', 'DIM_USER_0' in the NIfTI file.
    If the data size is inconsistent then it will only separate out the
    PVM_EncNReceivers into a 'DIM_COIL' dimension.

    Providing nfids will override the header derived shape of the data.

    :param path: Path to rawdata.job0
    :type path: pathlib.Path or str
    :param out: Output name and path
    :type out: pathlib.Path or str
    :param nfids: Number of FIDs expected in file, defaults to None
    :type nfids: int, optional
    """
    if isinstance(path, str):
        path = Path(path)
    if isinstance(out, str):
        out = Path(out)

    # Read parameters
    # From Method file
    method_info = _read_header_file_info(
        path.parent / 'method',
        ['PVM_EncNReceivers',
         'PVM_NAverages',
         'PVM_NRepetitions'],
        ['PVM_SpecSWH',
         'PVM_FrqRef',
         'PVM_SpecMatrix',
         'PVM_FrqWorkOffsetPpm'])

    # From visu_pars file
    visu_pars_info = _read_header_file_info(
        path.parent / 'visu_pars',
        [],
        ['VisuCoreSize'])

    # Check that there is a consistent FID length reported
    if visu_pars_info['VisuCoreSize'] != method_info['PVM_SpecMatrix']:
        print(
            'Inconsistent header information. '
            f'VisuCoreSize ({visu_pars_info["VisuCoreSize"]}) and '
            f'PVM_SpecMatrix ({method_info["PVM_SpecMatrix"]}) do not match. '
            'Using VisuCoreSize as FID length.')
        npoints = int(visu_pars_info['VisuCoreSize'])
    else:
        npoints = int(method_info['PVM_SpecMatrix'])

    # Load data using either the passed number of FIDs or header size
    if nfids is None:
        data = _load_array(path, npoints=npoints)
    else:
        data = _load_array(path, nfids=nfids)
        npoints = np.size(data,1)
        
    # Correct PVM_FrqWorkOffsetPpm if requested
    if correct_work_offset\
            and np.abs(method_info['PVM_FrqWorkOffsetPpm']) > 0:
        offset_hz = method_info['PVM_FrqRef']\
              * method_info['PVM_FrqWorkOffsetPpm']
        data = _correct_offset(
            data,
            1 / method_info['PVM_SpecSWH'],
            offset_hz)

    # reshape
    if nfids is None:
        proposed_shape = (
            method_info['PVM_NRepetitions'],
            method_info['PVM_NAverages'],
            method_info['PVM_EncNReceivers'],
            npoints)
    else:
         proposed_shape = (
            nfids,
            npoints)   

    if np.prod(proposed_shape) == np.prod(data.shape):
        if nfids is None:
            data_shaped = data.reshape(
                (method_info['PVM_NRepetitions'],
                 method_info['PVM_NAverages'],
                 method_info['PVM_EncNReceivers'],
                 npoints))
            tags = ['DIM_COIL', 'DIM_DYN', 'DIM_USER_0']
        else:
            data_shaped = data.reshape(
                (-1,
                 nfids,
                 npoints))
            tags = ['DIM_DYN', 'DIM_USER_0']        
    else:
        data_shaped = data.reshape(
            (-1,
             method_info['PVM_EncNReceivers'],
             npoints))
        tags = ['DIM_COIL', 'DIM_USER_0']

    data_shaped = data_shaped.T.reshape((1, 1, 1,) + data_shaped.T.shape)

    # Check output location exists
    if not out.parent.exists():
        out.parent.mkdir(parents=True)

    create_nmrs.gen_nifti_mrs(
        data_shaped,
        1 / method_info['PVM_SpecSWH'],
        method_info['PVM_FrqRef'],
        dim_tags=tags
    ).save(out)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert a Bruker rawdata.job0 file to NIfTI-MRS')
    parser.add_argument(
        'file',
        type=Path,
        help='path to rawdata.job0 file')
    parser.add_argument(
        'output',
        type=Path,
        help='name and path of output')
    parser.add_argument(
        '--nfid',
        required=False,
        type=int,
        help='Specified number of FIDs')
    parser.add_argument(
        '--no-work-offset-correction',
        action="store_false",
        dest="correct_work_offset",
        help='Disable the frequency shift correction '
             '(PVM_FrqWorkOffsetPpm).')
    args = parser.parse_args()

    main(
        args.file,
        args.output,
        nfids=args.nfid,
        correct_work_offset=args.correct_work_offset)
