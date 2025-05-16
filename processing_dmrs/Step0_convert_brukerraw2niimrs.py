#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert bruker raw data job0 file to fsl mrs

Last changed May 2025
@author: Malte O
"""
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils import preproc as proc
from fsl_mrs.core import nifti_mrs as ntools
from fsl_mrs.utils import preproc
from fsl_mrs.core import NIFTI_MRS
from scipy.io import loadmat
from custom_functions import *
from mrs_plots import *
import importlib
from bids_structure import *

import sys
sys.path.append('../common/nifti_mrs_from_raw')
sys.path.append('../')
import importlib

import fsl_mrs.utils.mrs_io as mrs_io
from fsl_mrs.utils.preproc import nifti_mrs_proc as proc


from fsl_mrs.core import nifti_mrs as ntools

from nifti_mrs_from_raw_pv3 import main as bruker2niimrs
import nifti_mrs_from_raw_pv3 as bruker2fsl

from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['mrs_plots'])

def coil_combine_bruker_header(data, cfg, return_uncombined_data =False):
    BRUKERparamnames = ["PVM_ArrayPhase",
                        "PVM_EncChanScaling"];

    path_to_data = cfg['data_path']
    water_reference_seqeunce_number = cfg['water_ref_seq_number']

    hdr = _read_header_file_info(f'{path_to_data}/{water_reference_seqeunce_number}/method', [], BRUKERparamnames)

    rxarrayphases = hdr[BRUKERparamnames[0]]
    scalingloops = hdr[BRUKERparamnames[1]]

    coil_factors = scalingloops * np.exp(1j * 2 * np.pi * rxarrayphases / 360)

    data_factored = data.copy()
    data_cc = data.copy(remove_dim='DIM_COIL')
    coil_dim = data.dim_position('DIM_COIL')

    if coil_dim==4 and len(data.shape)==5:
        data_factored[:] = data[:, :, :, :, :] * coil_factors[None, None, None, None, :]
    elif coil_dim == 4 and len(data.shape) == 6:
        data_factored[:] = data[:, :, :, :, :,:] * coil_factors[None, None, None, None, :, None]
    elif coil_dim==5 and len(data.shape) == 6:
        data_factored[:] = data[:, :, :, :, :, :] * coil_factors[None, None, None, None, :, None]
    else:
        return 'Unimplemented coil dimension.'
    if return_uncombined_data:
        return data_factored

    data_cc[:] = data_factored[:].sum(axis=coil_dim)
    return data_cc

def Step0_convert_bruker(subj_list, cfg):
    path_to_data = cfg['data_path']
    water_reference_seqeunce_number = cfg['water_ref_seq_number']
    min_seq_number = cfg['min_seq_number']
    max_seq_number = cfg['max_seq_number']
    fixed_phase_shift = cfg['fixed_phase_shift']

    ppmlim_outlier = (0.2, 4.3)
    subject_name  = subj_list[0]
    if len(subj_list) > 1:
        return 'More than one subject not yet implemented! Aborting.'

    # PROCESS WATER REFERENCE DATA
    ## Convert bruker raw data to nifti mrs
    bruker2niimrs(f'{path_to_data}/{water_reference_seqeunce_number}/rawdata.job0',
                  f'{path_to_data}/{water_reference_seqeunce_number}/rawdata_nifti')

    ref_data = mrs_io.read_FID(f'{path_to_data}/{water_reference_seqeunce_number}/rawdata_nifti.nii.gz')

    ## crop water reference along trivial dimensions
    if ref_data.shape[-1] == 1:
        ref_data = ntools.reshape(ref_data, (ref_data.shape[-3], ref_data.shape[-1], ref_data.shape[-2]))
    for i, dim in enumerate(ref_data.dim_tags):
        if ref_data.shape[4 + i] == 1:
            ref_data_crop = ref_data.copy(remove_dim=4 + i)
    ref_data = ref_data_crop

    ## for SPECIAL sequence: split even and odd repetitions
    even_reps = [2 * n for n in range(ref_data.shape[5] // 2)]
    ref_data_even, ref_data_odd = ntools.split(ref_data, 5, even_reps)

    avg_ref_data_even = proc.average(ref_data_even, 'DIM_USER_0')
    avg_ref_data_odd = proc.average(ref_data_odd, 'DIM_USER_0')

    ## coil combination
    if cfg['coil_combination_method'] == 'FSL MRS':
        avg_ref_data_even_cc = proc.coilcombine(avg_ref_data_even, reference=avg_ref_data_even)
        avg_ref_data_odd_cc = proc.coilcombine(avg_ref_data_odd, reference=avg_ref_data_odd)
    elif cfg['coil_combination_method'] == 'Bruker header info':
        avg_ref_data_even_cc = coil_combine_bruker_header(avg_ref_data_even,cfg)
        avg_ref_data_odd_cc = coil_combine_bruker_header(avg_ref_data_odd,cfg)
    else:
        return 'Please choose valid coil combination method. Currently available: \'FSL MRS\' or \'Bruker header info\'. Aborting.'

    ## fixed phase shift
    avg_ref_data_even_cc_shift = proc.apply_fixed_phase(avg_ref_data_even_cc, fixed_phase_shift)
    avg_ref_data_odd_cc_shift = proc.apply_fixed_phase(avg_ref_data_odd_cc, fixed_phase_shift)

    # SPECIAL: add even and odd averages
    avg_ref_data_combined_cc_shift = proc.add(avg_ref_data_even_cc_shift, avg_ref_data_odd_cc_shift)

    # PROCESS METABOLITE DATA
    ## convert raw data and load
    data = []
    for seq_number in range(min_seq_number, max_seq_number + 1):
        bruker2niimrs(f'{path_to_data}/{seq_number}/rawdata.job0', f'{path_to_data}/{seq_number}/rawdata_nifti')
        data.append(mrs_io.read_FID(f'{path_to_data}/{seq_number}/rawdata_nifti.nii.gz'))

    ## remove trivial dimension
    for j, this_data in enumerate(data):
        if this_data.shape[-1] == 1:
            this_data = ntools.reshape(this_data, (this_data.shape[-3], this_data.shape[-1], this_data.shape[-2]))
        for i, dim in enumerate(this_data.dim_tags):
            if this_data.shape[4 + i] == 1:
                data[j] = this_data.copy(remove_dim=dim)

    ## for SPECIAL processing: split even and odd repetitions
    data_even = []
    data_odd = []
    for this_data in data:
        even_reps = [2 * n for n in range(this_data.shape[5] // 2)]
        this_data_even, this_data_odd = ntools.split(this_data, 5, even_reps)
        data_even.append(this_data_even)
        data_odd.append(this_data_odd)

    ## coil combination
    data_even_cc = []
    data_odd_cc = []

    if cfg['coil_combination_method'] == 'FSL MRS':
        for this_data in data_even:
            data_even_cc.append(proc.coilcombine(this_data, reference=avg_ref_data_even, figure=False))
        for this_data in data_odd:
            data_odd_cc.append(proc.coilcombine(this_data, reference=avg_ref_data_odd, figure=False))
    elif cfg['coil_combination_method'] == 'Bruker header info':
        for this_data in data_even:
            data_even_cc.append(coil_combine_bruker_header(this_data, cfg))
        for this_data in data_odd:
            data_odd_cc.append(coil_combine_bruker_header(this_data,cfg ))
    else:
        return 'Please choose valid coil combination method. Currently available: \'FSL MRS\' or \'Bruker header info\'. Aborting.'

    ## shift by fixed shit
    data_even_cc_shift = []
    data_odd_cc_shift = []

    for this_data in data_even_cc:
        data_even_cc_shift.append(proc.apply_fixed_phase(this_data, fixed_phase_shift))

    for this_data in data_odd_cc:
        data_odd_cc_shift.append(proc.apply_fixed_phase(this_data, fixed_phase_shift))

    ## align averages
    data_even_cc_shift_align = []
    data_odd_cc_shift_align = []

    for this_data in data_even_cc_shift:
        data_even_cc_shift_align.append(proc.align(this_data, 'DIM_USER_0', ppmlim=(0.2, 4.2), figure=False))

    for this_data in data_odd_cc_shift:
        data_odd_cc_shift_align.append(proc.align(this_data, 'DIM_USER_0', ppmlim=(0.2, 4.2), figure=False))

    ## outlier removal
    data_even_cc_shift_align_unlike_good_ind = []
    data_odd_cc_shift_align_unlike_good_ind = []

    data_even_cc_shift_align_unlike_bad_ind = []
    data_odd_cc_shift_align_unlike_bad_ind = []

    for this_data in data_even_cc_shift_align:
        this_unlike = remove_unlike(this_data, ppmlim=ppmlim_outlier, return_indices=True)
        data_even_cc_shift_align_unlike_good_ind.append(this_unlike[0])
        data_even_cc_shift_align_unlike_bad_ind.append(this_unlike[1])
    for this_data in data_odd_cc_shift_align:
        this_unlike = remove_unlike(this_data, ppmlim=ppmlim_outlier, return_indices=True)
        data_odd_cc_shift_align_unlike_good_ind.append(this_unlike[0])
        data_odd_cc_shift_align_unlike_bad_ind.append(this_unlike[1])

    data_combined_cc_shift_align_unlike = []
    for j, this_data in enumerate(data_even_cc_shift_align):
        these_good_data = []
        for i in range(this_data.shape[-1]):
            # only consider data where both even and odd version is not an outlier
            if (i in data_even_cc_shift_align_unlike_good_ind[j]) and (i in data_odd_cc_shift_align_unlike_good_ind[j]):
                # take out data to add
                scrap, this_even_rep = ntools.split(data_even_cc_shift_align[j], 4, [i])
                scrap, this_odd_rep = ntools.split(data_odd_cc_shift_align[j], 4, [i])
                # add data and append
                these_good_data.append(proc.add(this_even_rep, this_odd_rep))
        # stack along repetition dimension
        data_combined_cc_shift_align_unlike.append(ntools.merge(these_good_data, 'DIM_USER_0'))
    print('Warning: We are adding two averages and not dividing by 2, maybe confusing quantification.')

    ## averaging data
    data_combined_cc_shift_align_unlike_av = []

    for this_data in data_combined_cc_shift_align_unlike:
        data_combined_cc_shift_align_unlike_av.append(proc.average(this_data, 'DIM_USER_0'))

    ## eddy current correction
    data_combined_cc_shift_align_unlike_av_ecc = []

    for this_data in data_combined_cc_shift_align_unlike_av:
        data_combined_cc_shift_align_unlike_av_ecc.append(proc.ecc(this_data, avg_ref_data_combined_cc_shift))

    print('Truncation etc needed? Double check Jessies script.')

    limits = [-0.15, 0.15]
    limunits = 'ppm'

    data_combined_cc_shift_align_unlike_av_ecc_watersupp = [proc.remove_peaks(this_data, limits, limit_units=limunits)
                                                            for this_data in data_combined_cc_shift_align_unlike_av_ecc]

    data_combined_cc_shift_align_unlike_av_ecc_watersupp_refshift = [
        proc.shift_to_reference(this_data, 3.027, (2.9, 3.1))
        for this_data in data_combined_cc_shift_align_unlike_av_ecc_watersupp]

    data_combined_cc_shift_align_unlike_av_ecc_watersupp_refshift_phased = [proc.phase_correct(this_data, (2.9, 3.1))
                                                                            for this_data in
                                                                            data_combined_cc_shift_align_unlike_av_ecc_watersupp_refshift]
    avg_ref_data_combined_cc_shift_phased = proc.phase_correct(avg_ref_data_combined_cc_shift, (4.55, 4.7), hlsvd=False)

    # add meta data
    doc = dict()
    doc['Bvalue'] = 'B value in ms / mu m^2.'
    doc['Delta'] = 'Diffusion gratient duration in ms.'
    doc['DiffusionTime'] = 'Diffusion time in ms.'
    doc['MixingTime'] = 'Mixing time in ms.'

    ## add water reference meta data
    meta_data = _read_header_file_info(f'{path_to_data}/{water_reference_seqeunce_number}/method',
                                       ['Bvalue', 'Delta', 'DiffusionTime', 'MixingTime'], [])
    for key in meta_data:
        avg_ref_data_combined_cc_shift_phased.add_hdr_field(key, meta_data[key], doc=doc[key])

    ## add metabolite meta data
    for seq_number in range(min_seq_number, max_seq_number + 1):
        meta_data = _read_header_file_info(f'{path_to_data}/{seq_number}/method',
                                           ['Bvalue', 'Delta', 'DiffusionTime', 'MixingTime'], [])
        for key in meta_data:
            data_combined_cc_shift_align_unlike_av_ecc_watersupp_refshift_phased[
                seq_number - min_seq_number].add_hdr_field(key, meta_data[key], doc=doc[key])

    # SAVE DATA

    ## create bids structure
    bids_strc = create_bids_structure(subj=subject_name, sess=1, datatype='dmrs', root=path_to_data,
                                      folderlevel='derivatives', workingdir=cfg['prep_foldername'] )

    if not os.path.exists(bids_strc.get_path()):
        os.makedirs(bids_strc.get_path())

    for i, this_data in enumerate(data_combined_cc_shift_align_unlike_av_ecc_watersupp_refshift_phased):
        this_data.save(bids_strc.get_path(f'seq_{min_seq_number + i}_dmrs.nii.gz'))
    avg_ref_data_combined_cc_shift_phased.save(
        bids_strc.get_path(f'seq_{water_reference_seqeunce_number}_ref_dmrs.nii.gz'))
    return "All data preprocessed successfully."

def Step0_combine_diffusion_times(subj_list, cfg):
    for subject_name in subj_list:
        for ses_no in [1]:
            bids_strc = create_bids_structure(subj=subject_name, sess=ses_no, datatype='dmrs', root=cfg['data_path'],
                                              folderlevel='derivatives', workingdir=cfg['prep_foldername'])

            print('Warning: Some manual removal of sequences when combining them.')
            niftis = [mrs_io.read_FID(f'{bids_strc.get_path()}/{subject_name}_ses-0{ses_no}_seq_{seq_no}_dmrs.nii.gz')
                      for seq_no in [cfg['min_seq_number']] + [i for i in range(cfg['min_seq_number']+2, cfg['max_seq_number'] + 1)]]

            print('Found diffusion times ', np.unique([nifti.hdr_ext['DiffusionTime']['Value'] for nifti in niftis]))
            if cfg['diffusion_times'] == 'all':
                DiffusionTimes = np.unique([nifti.hdr_ext['DiffusionTime']['Value'] for nifti in niftis])
            else:
                DiffusionTimes = cfg['diffusion_times']

            for DiffusionTime in DiffusionTimes:
                these_niftis = []
                bvals = []
                for nifti in niftis:
                    if nifti.hdr_ext['DiffusionTime']['Value'] == DiffusionTime:
                        nifti.set_dim_tag(4, 'DIM_USER_0')
                        these_niftis.append(nifti)
                        bvals.append(nifti.hdr_ext['Bvalue']['Value'])
                nifti_dict = dict(zip(bvals, these_niftis))
                these_niftis_sorted = [nifti_dict[bval] for bval in sorted(bvals)]
                big_nifti = ntools.merge(these_niftis_sorted, dimension=4)
                big_nifti.set_dim_tag(
                    'DIM_USER_0',
                    'DIM_USER_0',
                    'b-value',
                    {'b_value': {'Value': list(1e-3*np.array(sorted(bvals))), 'Description': 'b-value in ms.μm^-2'}})
                big_nifti.save(
                    f'{bids_strc.get_path()}/{subject_name}_ses-0{ses_no}_TD_{int(DiffusionTime)}_dmrs.nii.gz')
    return 'All niftis grouped by diffusion times'

            # from https://git.fmrib.ox.ac.uk/fsl/fsl_mrs/-/blob/master/fsl_mrs/utils/preproc/nifti_mrs_proc.py

def remove_unlike(data, ppmlim=None, sdlimit=1.96, niter=2, figure=False, report=None, return_indices=False):
    '''Remove unlike dynamics operating on DIM_DYN

    :param NIFTI_MRS data: Data to truncate or pad
    :param figure: True to show figure.
    :param report: Provide output location as path to generate report

    :return: Data passing likeness criteria.
    :return: Data failing likness criteria
    '''
    if data.shape[:3] != (1, 1, 1):
        raise OnlySVS("remove_unlike only specified for SVS data")

    if data.ndim > 5:
        raise ValueError('remove_unlike only makes sense for a single dynamic dimension. Combined coils etc. first')
    elif data.ndim < 5:
        raise ValueError('remove_unlike only makes sense for data with a dynamic dimension')

    goodFIDs, badFIDs, gIndicies, bIndicies, metric = \
        preproc.identifyUnlikeFIDs(data[0, 0, 0, :, :].T,
                                   data.bandwidth,
                                   data.spectrometer_frequency[0],
                                   nucleus=data.nucleus[0],
                                   ppmlim=ppmlim,
                                   sdlimit=sdlimit,
                                   iterations=niter,
                                   shift=True)

    if figure or report:
        from fsl_mrs.utils.preproc.unlike import identifyUnlikeFIDs_report
        fig = identifyUnlikeFIDs_report(goodFIDs,
                                        badFIDs,
                                        gIndicies,
                                        bIndicies,
                                        metric,
                                        data.bandwidth,
                                        data.spectrometer_frequency[0],
                                        nucleus=data.nucleus[0],
                                        ppmlim=ppmlim,
                                        sdlimit=sdlimit,
                                        html=report)
        if figure:
            fig.show()

    if return_indices:
        return gIndicies, bIndicies

    # goodFIDs = np.asarray(goodFIDs).T
    # goodFIDs = goodFIDs.reshape([1, 1, 1] + list(goodFIDs.shape))

    if len(badFIDs) > 0:
        bad_out, good_out = ntools.split(
            data,
            data.dim_tags[0],
            gIndicies)
    else:
        good_out = data.copy()

    good_out.add_hdr_field(
        f'{data.dim_tags[0]} Indices',
        gIndicies,
        doc=f"Data's original index values in the {data.dim_tags[0]} dimension")

    if len(badFIDs) > 0:
        bad_out.add_hdr_field(
            f'{data.dim_tags[0]} Indices',
            bIndicies,
            doc=f"Data's original index values in the {data.dim_tags[0]} dimension")
    else:
        bad_out = None

    # Update processing prov
    processing_info = f'{__name__}.remove_unlike, '
    if ppmlim is None:
        processing_info += 'ppmlim=None, '
    else:
        processing_info += f'ppmlim={ppmlim}, '
    processing_info += f'sdlimit={sdlimit}, '
    processing_info += f'niter={niter}.'

    update_processing_prov(good_out, 'Outlier removal', processing_info)

    return good_out, bad_out


from fsl_mrs import __version__


def update_processing_prov(nmrs_obj: NIFTI_MRS, method, details):
    """Insert appropriate processing provenance information into the
    NIfTI-MRS header extension.

    :param nmrs_obj: NIFTI-MRS object which has been modified
    :type nmrs_obj: fsl_mrs.core.NIFTI_MRS
    :param method: [description]
    :type method: str
    :param details: [description]
    :type details: str
    """
    # 1. Check for ProcessingApplied key and create if not present
    if 'ProcessingApplied' in nmrs_obj.hdr_ext:
        current_processing = nmrs_obj.hdr_ext['ProcessingApplied']
    else:
        current_processing = []

    # 2. Form object to append.
    prov_dict = {
        'Time': datetime.now().isoformat(sep='T', timespec='milliseconds'),
        'Program': 'FSL-MRS',
        'Version': __version__,
        'Method': method,
        'Details': details}

    # 3. Append
    current_processing.append(prov_dict)
    nmrs_obj.add_hdr_field('ProcessingApplied', current_processing)


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
                method_values[match[1]] = np.array(
                    methodlines[idx+1].split(' ')).astype('float')

    return method_values
