import os
import copy
import sys
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import rcParams
import shared_variables

# Set font to Arial and adjust font sizes
rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 14,  # General font size
    'axes.titlesize': 18,  # Title font size
    'axes.labelsize': 16,  # Axis label font size
    'xtick.labelsize': 14,  # X-axis tick label size
    'ytick.labelsize': 14,  # Y-axis tick label size
    'legend.fontsize': 14  # Legend font size
})

# The path to the LOGS directory
LOG_DIR = '/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/LOGS'

# List the directories to the resource logging files for each sample
# SAMPLE_LIST = ["sample_1000", "sample_2000", "sample_3000", "sample_4000", "sample_5000"]

def parse_wall_clock_time(line):
    # Extract the time part after the last mention of 'time'
    time_part = line.split("):")[-1].strip()
    
    # Split the time part by colons to get hours, minutes, and seconds if present
    time_parts = time_part.split(":")
    
    # Initialize hours, minutes, seconds to 0
    hours, minutes, seconds = 0, 0, 0
    
    # Clean up and parse each part
    if len(time_parts) == 3:  # h:mm:ss or h:mm:ss.ss
        hours = float(re.sub(r'[^\d.]', '', time_parts[0]))  # Remove non-numeric characters
        minutes = float(re.sub(r'[^\d.]', '', time_parts[1]))
        seconds = float(re.sub(r'[^\d.]', '', time_parts[2]))

    elif len(time_parts) == 2:  # m:ss or m:ss.ss
        minutes = float(re.sub(r'[^\d.]', '', time_parts[0]))
        seconds = float(re.sub(r'[^\d.]', '', time_parts[1]))

    # Calculate total time in seconds
    total_seconds = seconds + (minutes * 60) + (hours * 3600)
    return total_seconds

def plot_metric_by_step_adjusted(sample_resource_dict, metric, ylabel, title, filename, divide_by_cpu=False):
    """
    Plots the specified metric for each sample and step.
    
    Args:
    - sample_resource_dict: Dictionary containing sample resource data.
    - metric: The metric to plot (e.g., 'user_time', 'system_time').
    - ylabel: Label for the y-axis.
    - title: Title of the plot.
    - filename: Path to save the plot.
    - divide_by_cpu: Boolean indicating if the metric should be divided by 'percent_cpu'.
    """
    # Extract the metric data for each sample and step
    samples = sorted(sample_resource_dict.keys())
    steps = []
    for sample in samples:
        # print(sample_resource_dict[sample].keys())
        for key in sample_resource_dict[sample].keys():
            steps.append(key)
    
    # print(steps)
    # print(sample_resource_dict)
    
    metric_data = {sample: [
        (sample_resource_dict[sample][step][metric] / sample_resource_dict[sample][step]['percent_cpu'] 
         if divide_by_cpu else sample_resource_dict[sample][step][metric])
        for step in steps
    ] for sample in samples}

    # Plotting the figure with each step on the x-axis and the specified metric on the y-axis
    x = np.arange(len(steps))  # x locations for the groups
    bar_width = 0.02  # Width of the bars
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot each sample as a bar group, with each sample's bars slightly offset
    for i, sample in enumerate(samples):
        ax.bar(x + i * bar_width, metric_data[sample], bar_width, label=f'{sample.split("_")[0]} Cells')

    # Labeling the plot
    ax.set_xlabel('Step')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x + bar_width * (len(samples) - 1) / 2)
    ax.set_xticklabels(steps, rotation=45)
    plt.figlegend(loc='center right', bbox_to_anchor=(0.98, 0.5), fontsize=8, ncol=1)

    # Show plot
    plt.tight_layout(rect=[0, 0.06, 0.85, 0.90])
    plt.savefig(filename, dpi=200)
    print(f'Saved {filename.split("/")[-1]}')

def plot_total_metric_by_sample(sample_resource_dict, metric, ylabel, title, filename, divide_by_cpu=False):
    """
    Plots the total specified metric across all steps for each sample.
    
    Args:
    - sample_resource_dict: Dictionary containing sample resource data.
    - metric: The metric to plot (e.g., 'user_time', 'system_time').
    - ylabel: Label for the y-axis.
    - title: Title of the plot.
    - filename: Path to save the plot.
    - divide_by_cpu: Boolean indicating if the metric should be divided by 'percent_cpu'.
    """
    # Calculate the total metric for each sample across all steps
    samples = sorted(sample_resource_dict.keys())
    

    
    if metric == 'max_ram':
        total_metric_data = [max(sample_resource_dict[sample][step][metric] for step in sample_resource_dict[sample]) for sample in samples]
    if metric == 'percent_cpu':
        total_metric_data = [sum(sample_resource_dict[sample][step][metric] for step in sample_resource_dict[sample]) / len(sample_resource_dict[sample]) for sample in samples]
    else:
        total_metric_data = []
        for sample in samples:
            total = sum(
                sample_resource_dict[sample][step][metric] / sample_resource_dict[sample][step]['percent_cpu'] 
                if divide_by_cpu else sample_resource_dict[sample][step][metric]
                for step in sample_resource_dict[sample]
            )
            total_metric_data.append(total)

    # Plotting the total metric for each sample
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(samples))  # x locations for the samples
    bar_width = 0.15  # Width of the bars


    ax.bar(x, total_metric_data, bar_width, color='blue', label=metric)

    # Labeling the plot
    ax.set_xlabel('Number of Cells')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels([sample.split("_")[0] for sample in samples], rotation=45)

    # Show plot
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    print(f'Saved {filename.split("/")[-1]}')

if __name__ == '__main__':
    # Define a dictionary to hold the sample names with their resource requirements for each step in the pipeline
    sample_resource_dict = {}
    
    samples = ["SCENIC_PLUS_LOGS"]

    # List of files in the LOGS folder matching the sample names
    sample_list = [
        sample_dir for sample_dir in os.listdir(LOG_DIR)
        if any(rep in sample_dir for rep in samples)
    ]
    print(f'\nSample list {sample_list}')

    # For each  
    for sample_log_dir in os.listdir(LOG_DIR):
        if sample_log_dir in samples:
            print(f'\nAnalyzing {sample_log_dir}')
            
            sample_resource_dict[sample_log_dir] = {}
            
            # Find each step log file for the sample
            for file in os.listdir(f'{LOG_DIR}/{sample_log_dir}'):
                if file.endswith(".log") and "step" in file.lower():
                    print(f'\t{file}')
                    pipeline_step = file.split(".")[0]
                    sample_resource_dict[sample_log_dir][pipeline_step] = {
                        "user_time": 0,
                        "system_time": 0,
                        "percent_cpu": 0,
                        "wall_clock_time": 0,
                        "max_ram": 0
                    }

                    # Extract each relevant resource statistic for the sample step and save it in a dictionary
                    with open(f'{LOG_DIR}/{sample_log_dir}/{file}', 'r') as log_file:
                        for line in log_file:
                            if 'User time' in line:
                                sample_resource_dict[sample_log_dir][pipeline_step]["user_time"] = float(line.split(":")[-1])
                            if 'System time' in line:
                                sample_resource_dict[sample_log_dir][pipeline_step]["system_time"] = float(line.split(":")[-1])
                            if 'Percent of CPU' in line:
                                sample_resource_dict[sample_log_dir][pipeline_step]["percent_cpu"] = float(line.split(":")[-1].split("%")[-2])
                            if 'wall clock' in line:
                                sample_resource_dict[sample_log_dir][pipeline_step]["wall_clock_time"] = parse_wall_clock_time(line)
                            if 'Maximum resident set size' in line:
                                kb_per_gb = 1048576
                                sample_resource_dict[sample_log_dir][pipeline_step]["max_ram"] = (float(line.split(":")[-1]) / kb_per_gb)

    # Plot the resource requirements by step for each sample
    plot_metric_by_step_adjusted(
        sample_resource_dict=sample_resource_dict,
        metric='user_time',
        ylabel='User Time (s) / Percent CPU Usage',
        title='User Time / Percent CPU Usage by Step for Each Sample',
        filename=f'{shared_variables.results_dir}/resource_analysis/Step_User_Time_Summary.png',
        divide_by_cpu=True
    )

    plot_metric_by_step_adjusted(
        sample_resource_dict=sample_resource_dict,
        metric='system_time',
        ylabel='System Time (s) / Percent CPU Usage',
        title='System Time / Percent CPU Usage by Step for Each Sample',
        filename=f'{shared_variables.results_dir}/resource_analysis/Step_System_Time.png',
        divide_by_cpu=True
    )

    plot_metric_by_step_adjusted(
        sample_resource_dict=sample_resource_dict,
        metric='wall_clock_time',
        ylabel='Wall Clock Time (s)',
        title='Wall Clock Time by Step for Each Sample',
        filename=f'{shared_variables.results_dir}/resource_analysis/Step_Wall_Clock_Time.png',
        divide_by_cpu=False
    )

    plot_metric_by_step_adjusted(
        sample_resource_dict=sample_resource_dict,
        metric='max_ram',
        ylabel='Max RAM Usage (GB)',
        title='Max RAM usage by Step for Each Sample',
        filename=f'{shared_variables.results_dir}/resource_analysis/Step_Max_Ram.png',
        divide_by_cpu=False
    )

    plot_metric_by_step_adjusted(
        sample_resource_dict=sample_resource_dict,
        metric='percent_cpu',
        ylabel='Percent CPU',
        title='Percent of the CPU Used',
        filename=f'{shared_variables.results_dir}/resource_analysis/Step_Percent_Cpu.png',
        divide_by_cpu=False
    )

    # Plot the resource requirements for running the entire pipeline
    plot_total_metric_by_sample(
        sample_resource_dict=sample_resource_dict,
        metric='user_time',
        ylabel='Total User Time / Percent CPU Usage',
        title='Total User Time / Percent CPU Usage for Each Sample',
        filename=f'{shared_variables.results_dir}/resource_analysis/Total_User_Time.png',
        divide_by_cpu=True
    )

    plot_total_metric_by_sample(
        sample_resource_dict=sample_resource_dict,
        metric='system_time',
        ylabel='Total System Time / Percent CPU Usage',
        title='Total System Time / Percent CPU Usage for Each Sample',
        filename=f'{shared_variables.results_dir}/resource_analysis/Total_System_Time.png',
        divide_by_cpu=True
    )

    plot_total_metric_by_sample(
        sample_resource_dict=sample_resource_dict,
        metric='wall_clock_time',
        ylabel='Wall Clock Time (s)',
        title='Total Wall Clock Time',
        filename=f'{shared_variables.results_dir}/resource_analysis/Total_Wall_Clock_Time.png',
        divide_by_cpu=False
    )

    plot_total_metric_by_sample(
        sample_resource_dict=sample_resource_dict,
        metric='max_ram',
        ylabel='Max RAM Usage (GB)',
        title='Max RAM usage',
        filename=f'{shared_variables.results_dir}/resource_analysis/Total_Max_Ram.png',
        divide_by_cpu=False
    )

    plot_total_metric_by_sample(
        sample_resource_dict=sample_resource_dict,
        metric='max_ram',
        ylabel='Max RAM Usage (GB)',
        title='Max RAM usage',
        filename=f'{shared_variables.results_dir}/resource_analysis/Total_Max_Ram.png',
        divide_by_cpu=False
    )

    plot_total_metric_by_sample(
        sample_resource_dict=sample_resource_dict,
        metric='percent_cpu',
        ylabel='Percent CPU',
        title='Average Percent of the CPU Used',
        filename=f'{shared_variables.results_dir}/resource_analysis/Total_Percent_Cpu.png',
        divide_by_cpu=False
    )

    summary_dict = {}

    for sample, step_dict in sample_resource_dict.items():
        print(sample)
        if sample not in summary_dict:
            summary_dict[sample] = {
                    "user_time": 0,
                    "system_time": 0,
                    "percent_cpu": [],
                    "wall_clock_time": 0,
                    "max_ram": []
                }
        for step, resource_dict in step_dict.items():
            print(f'\t{step}')
            for resource_name, resource_value in resource_dict.items():
                print(f'\t\t{resource_name}: {resource_value}')
                if resource_name == "percent_cpu":
                    summary_dict[sample][resource_name].append(round(resource_value,2))
                elif resource_name == "max_ram":
                    summary_dict[sample][resource_name].append(round(resource_value,2))
                else:
                    summary_dict[sample][resource_name] += round(resource_value,2)
        summary_dict[sample]["max_ram"] = max(summary_dict[sample]["max_ram"])
        summary_dict[sample]["percent_cpu"] = round(sum(summary_dict[sample]["percent_cpu"]) / len(summary_dict[sample]["percent_cpu"]),2)
        
    summary_df = pd.DataFrame(summary_dict)
    summary_df = summary_df.reindex(sorted(summary_df.columns), axis=1)
    print(summary_df.head())

    summary_df.to_csv(f'{shared_variables.results_dir}/resource_analysis/Resource_Summary.tsv', sep='\t')