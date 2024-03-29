import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DateFormatter
import pandas as pd
from pathlib import Path

import utils.io_routines as io

def plot_min_infected(mi_table, file):
    file_path = Path(file)
    file_path.parent.mkdir(parents=True, exist_ok=True)  
    fig, ax = plt.subplots()
    plt.grid(True, alpha=0.5)
    plt.fill_between(pd.to_datetime(mi_table['date']), 
             mi_table['min_n_true'], 
             mi_table['smoothed_cases'], 
             alpha=0.7,
             color='darkblue',
             label="min. undetected cases")
    plt.fill_between(pd.to_datetime(mi_table['date']), 
             mi_table['smoothed_cases'], 
             0, 
             alpha=0.7,
             color='darkred',
             label="reported cases")
    plt.ylabel("min. number of new infections")

    #ax.legend()
    plt.legend(bbox_to_anchor =(0.75, 1.15), ncol = 2)

    # set plot size
    fig.set_size_inches(10, 5)
    # Format date in x axis (month and year)
    date_form = DateFormatter("%b %Y")
    ax.xaxis.set_major_formatter(date_form)
    ax.xaxis.set_major_locator(MonthLocator())
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
    #plt.show()

    # bbox_inches = tight to take also labels of axis into account (otherwise it's cut)
    plt.savefig(file, bbox_inches='tight', dpi=500)

    # close/refresh plot
    plt.clf()
    plt.cla()
    plt.close()


def plot_phi_estimates(phi_table, smoothed_phi_table, file):
    file_path = Path(file)
    file_path.parent.mkdir(parents=True, exist_ok=True)  
    plot_table = pd.DataFrame({'date': pd.to_datetime(phi_table['date']),
                               'phi': phi_table['phi'],
                               'size': phi_table['sampleSize']/max(phi_table['sampleSize'])*100
                               })

    fig, ax = plt.subplots()
    plt.rc('text.latex', preamble=r'\usepackage{textgreek}')
    plt.grid(True, alpha=0.5)
    plt.scatter('date', 'phi', 
                s='size',
                c='black', 
                alpha=0.5,
                data=plot_table)
    plt.plot(pd.to_datetime(smoothed_phi_table['date']), 
             smoothed_phi_table['smoothed_phi'], 
             linewidth=5,
             alpha=0.7,
             c='white')
    plt.plot(pd.to_datetime(smoothed_phi_table['date']), 
             smoothed_phi_table['smoothed_phi'], 
             linewidth=3,
             c='darkblue')
    plt.ylabel(r'incidence correlate $\phi$')

    # set plot size
    fig.set_size_inches(10, 5)

    # Format date in x axis (month and year)
    date_form = DateFormatter("%b %Y")
    ax.xaxis.set_major_formatter(date_form)
    ax.xaxis.set_major_locator(MonthLocator())
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
    
    plt.savefig(file, bbox_inches='tight', dpi=500)

    # close/refresh plot
    plt.clf()
    plt.cla()
    plt.close()