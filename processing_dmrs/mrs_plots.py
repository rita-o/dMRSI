import numpy as np
import matplotlib.pyplot as plt
import os
from fsl_mrs.utils import plotting as splot
import matplotlib.gridspec as gridspec

def plot_spectrum(mrs_list, res_list=None, time_var=None, ppmlim=(-5, 15), proj='real', output_folder=None):
    """
    Plot dynamic MRS data over time using matplotlib for each time step.

    Args:
        mrs_list: list of MRS objects
        res_list: list of Results objects (optional)
        time_var: list of time variables (e.g., bvals for dwMRS)
        ppmlim: tuple (low, high) in ppm for x-axis range
        proj: string ('real', 'imag', 'abs', 'angle') to project the complex data
        output_folder: directory to save the generated plots

    Returns:
        None (plots saved as images)
    """
    if ppmlim is None and res_list is None:
        ppmlim = mrs_list[0].default_ppm_range
    elif ppmlim is None:
        ppmlim = res_list[0].ppmlim

    n = len(mrs_list)
    if time_var is None:
        time_var = np.arange(n)
    else:
        time_var = np.asarray(time_var)

    # Projection functions for the complex data
    proj_funcs = {'real': np.real,
                  'imag': np.imag,
                  'angle': np.angle,
                  'abs': np.abs}
    
    colors = {'data': 'gray', 'pred': 'red', 'base': 'blue', 'resid': 'lightgray'}
    
    # Prepare data to plot
    xaxis = mrs_list[0].getAxes()
    ydata = {'data': [proj_funcs[proj](mrs_list[i].get_spec()) for i in range(n)]}
    
    if res_list is not None:
        ydata['pred'] = [proj_funcs[proj](res_list[i].pred_spec) for i in range(n)]
        ydata['base'] = [proj_funcs[proj](splot.FID2Spec(res_list[i].baseline)) for i in range(n)]
        ydata['resid'] = [proj_funcs[proj](splot.FID2Spec(res_list[i].residuals)) for i in range(n)]
    

    # Generate and save static plots
    for i in range(n):
        t_str = np.array2string(time_var[i], floatmode='fixed', precision=1)

        plt.figure(figsize=(8, 6))
        for name in ydata:
            plt.plot(xaxis, ydata[name][i], label=f"{name} - {t_str}", color=colors[name], linewidth=1)
        
        plt.title(f'MRS Data at Time {t_str}')
        plt.xlabel('Chemical Shift (ppm)')
        plt.ylabel('Amplitude')
        plt.xlim(ppmlim)
        
        # Automatically adjust y-axis range
        data = np.asarray(ydata['data']).flatten()
        minval, maxval = np.min(data), np.max(data)
        ymin = minval - np.abs(minval) / 2
        ymax = maxval + maxval / 30
        plt.ylim([ymin, ymax])

        # Add legend
        plt.legend()

        # Save the figure
        if output_folder:
            # Ensure the output folder exists
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            plot_filename = os.path.join(output_folder, f'MRS_time_{t_str}.png')
            plt.savefig(plot_filename)
            plt.close()
            print(f"Plots saved in {output_folder}.")

# Example usage:
# Assuming mrs_list and res_list are defined as in your original code


def plot_fit(mrs, res, out=None, baseline=True, proj='real'):
    """Primary plotting function for FSL-MRS fits

    :param mrs: mrs object which has been fitted
    :type mrs: fsl_mrs.core.mrs.MRS
    :param res: Fitting results object
    :type res: fsl_mrs.utils.results.FitRes
    :param out: output figure filename, defaults to None
    :type out: str, optional
    :param baseline: optionally plot baseline, defaults to True
    :type baseline: bool, optional
    :param proj: 'real', 'imag', 'abs', or 'angle', defaults to 'real'
    :type proj: str, optional
    """

    def axes_style(plt, ppmlim, label=None, xticks=None):
        plt.xlim(ppmlim)
        plt.gca().invert_xaxis()
        plt.xlabel(label)
        plt.gca().set_xticks(xticks)
        plt.minorticks_on()
        plt.grid(True, axis='x', which='major', color='k', linestyle='--', linewidth=.3)
        plt.grid(True, axis='x', which='minor', color='k', linestyle=':', linewidth=.3)

    def doPlot(data, c='b', linewidth=1, linestyle='-', xticks=None):
        plt.plot(mrs.getAxes(), data, color=c, linewidth=linewidth, linestyle=linestyle)
        axes_style(plt, res.ppmlim, label='Chemical shift (ppm)', xticks=xticks)

    # Prepare data for plotting
    data = splot.FID2Spec(mrs.FID)
    pred = splot.FID2Spec(res.pred)
    if baseline is not None:
        baseline = splot.FID2Spec(res.baseline)

    first, last = mrs.ppmlim_to_range(ppmlim=res.ppmlim, shift=True)

    # turn to real numbers
    data = splot.data_proj(data, proj)
    pred = splot.data_proj(pred, proj)
    if baseline is not None:
        baseline = splot.data_proj(baseline, proj)

    if first > last:
        first, last = last, first

    m = min(data[first:last].min(), pred[first:last].min())
    M = max(data[first:last].max(), pred[first:last].max())
    ylim = (m - np.abs(M) / 10, M + np.abs(M) / 10)

    # Create the figure
    plt.figure(figsize=(9, 10))

    # Subplots
    gs = gridspec.GridSpec(2, 1,
                           height_ratios=[1, 20])

    plt.subplot(gs[0])
    # Start by plotting error
    if mrs.nucleus == '1H':
        xticks = np.arange(res.ppmlim[0], res.ppmlim[1], .2)
    else:
        xticks = np.round(np.linspace(res.ppmlim[0], res.ppmlim[1], 10), decimals=1)
    plt.plot(mrs.getAxes(), splot.data_proj(data - pred, proj), c='k', linewidth=1, linestyle='-')
    axes_style(plt, res.ppmlim, xticks=xticks)
    plt.gca().set_xticklabels([])

    plt.subplot(gs[1])

    doPlot(data, c='k', linewidth=.5, xticks=xticks)
    doPlot(pred, c='#cc0000', linewidth=1, xticks=xticks)
    if baseline is not None:
        doPlot(baseline, c='k', linewidth=.5, xticks=xticks)

    # plot y=0
    doPlot(data * 0, c='k', linestyle=':', linewidth=1, xticks=xticks)

    plt.legend(['data', 'model fit'])

    plt.tight_layout()
    plt.ylim(ylim)

    if out is not None:
        plt.savefig(out)

    return plt.gcf()
