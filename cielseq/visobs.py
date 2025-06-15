#!/usr/bin/env python3

import json
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.coordinates import EarthLocation, AltAz, ICRS, get_body
from astropy import units as u
from astropy.time import Time
from astroplan import Observer, FixedTarget
import numpy as np
import argparse
import datetime

def parse_arguments():
    """
    Parses command-line arguments for the plot date, JSON file path, and optional RA/Dec.
    Automatically handles --help to show usage information.
    """
    parser = argparse.ArgumentParser(
        description="Plot celestial object altitudes and twilight periods for a given date.",
        formatter_class=argparse.RawTextHelpFormatter # Keeps formatting for multiline help
    )
    parser.add_argument(
        '--date',
        type=str,
        help='Date for the plot in ISO-MM-DD format (e.g., 2025-06-11).\nDefaults to the current date if not provided.',
        default=None
    )
    parser.add_argument(
        '--json-file',
        type=str,
        help='Path to the JSON file containing celestial object data.\nDefaults to "fields.json" if not provided.\nIgnored if --ra and --dec are provided.',
        default='fields.json'
    )
    parser.add_argument(
        '--ra',
        type=str,
        help='Right Ascension (e.g., "10h30m00s" or "157.5d"). If provided, --dec must also be provided. Ignores --json-file.'
    )
    parser.add_argument(
        '--dec',
        type=str,
        help='Declination (e.g., "+30d00m00s" or "30d"). If provided, --ra must also be provided. Ignores --json-file.'
    )
    args = parser.parse_args()

    plot_date = None
    if args.date:
        try:
            plot_date = datetime.datetime.strptime(args.date, '%Y-%m-%d').date()
        except ValueError:
            print("Error: Invalid date format for --date. Please use ISO-MM-DD (e.g., 2025-06-11).")
            plot_date = datetime.date.today()
    else:
        plot_date = datetime.date.today()

    # Handle RA/Dec arguments
    ra_dec_provided = args.ra is not None or args.dec is not None
    if ra_dec_provided and (args.ra is None or args.dec is None):
        parser.error("Both --ra and --dec must be provided if one is used.")

    return plot_date, args.json_file, args.ra, args.dec

def load_json_data(file_path='fields.json'):
    """
    Loads celestial object data from a specified JSON file.
    Args:
        file_path (str): The path to the JSON file. Defaults to 'fields.json'.
    Returns:
        dict: The loaded JSON data, or None if an error occurred.
    """
    try:
        with open(file_path, 'r') as f:
            json_data = json.load(f)
        return json_data
    except FileNotFoundError:
        print(f"Error: '{file_path}' not found. Please make sure the file exists or provide a valid path.")
        return None
    except json.JSONDecodeError:
        print(f"Error: Could not decode '{file_path}'. Please check if the file contains valid JSON.")
        return None

def setup_observer_and_times(plot_date):
    """
    Defines the observer's location and generates a series of times for calculations.

    Args:
        plot_date (datetime.date): The date for which to generate the time series.

    Returns:
        tuple: A tuple containing:
            - observer (astroplan.Observer): The observer object.
            - times_utc (astropy.time.Time): A series of UTC times for the 24-hour period.
            - times_local (astropy.time.Time): A series of local times for the 24-hour period.
            - time_start_utc (astropy.time.Time): The start time (UTC) of the plot.
            - time_end_utc (astropy.time.Time): The end time (UTC) of the plot.
            - time_offset_local (astropy.units.quantity.Quantity): The local time offset from UTC.
            - tomorrow_date (datetime.date): The date for the next day.
    """
    # Define the observer's location (Königstein im Taunus, Germany)
    location = EarthLocation(lat=50.1833 * u.deg, lon=8.4667 * u.deg, height=400 * u.m)
    observer = Observer(location=location, name="Königstein im Taunus")

    # Construct the start and end times for a 24-hour cycle (noon to noon UTC)
    time_start_str = f"{plot_date.year}-{plot_date.month:02d}-{plot_date.day:02d}T12:00:00"
    tomorrow_date = plot_date + datetime.timedelta(days=1)
    time_end_str = f"{tomorrow_date.year}-{tomorrow_date.month:02d}-{tomorrow_date.day:02d}T12:00:00"

    time_start_utc = Time(time_start_str, format='isot', scale='utc')
    time_end_utc = Time(time_end_str, format='isot', scale='utc')

    # Generate a series of times for calculations (every 5 minutes in UTC)
    delta_t = 5 * u.min
    times_utc = Time(np.arange(time_start_utc.jd, time_end_utc.jd, delta_t.to(u.day).value), format='jd', scale='utc')

    # Determine local time offset (CEST for summer, CET for winter)
    if plot_date.month >= 3 and plot_date.month <= 10: # Rough estimation for DST period in Northern Hemisphere
        time_offset_local = 2 * u.hour # Central European Summer Time (CEST)
    else:
        time_offset_local = 1 * u.hour # Central European Time (CET)

    times_local = times_utc + time_offset_local # Calculate corresponding local times

    return observer, times_utc, times_local, time_start_utc, time_end_utc, time_offset_local, tomorrow_date

def create_targets(json_data, times_utc, observer_location, custom_ra=None, custom_dec=None):
    """
    Creates FixedTarget objects from the loaded JSON data or from custom RA/Dec, and adds the Moon.

    Args:
        json_data (dict): Dictionary containing celestial object coordinates and names.
        times_utc (astropy.time.Time): A series of UTC times used for Moon's ephemeris.
        observer_location (astropy.coordinates.EarthLocation): The observer's location.
        custom_ra (str, optional): Right Ascension string for a custom target.
        custom_dec (str, optional): Declination string for a custom target.

    Returns:
        list: A list of astroplan.FixedTarget objects.
    """
    targets = []

    if custom_ra and custom_dec:
        try:
            coords = ICRS(ra=custom_ra, dec=custom_dec)
            # Directly set the name to the formatted RA/Dec string for the legend
            target_name = f"RA: {custom_ra}, Dec: {custom_dec}"
            target = FixedTarget(coord=coords, name=target_name)
            targets.append(target)
        except ValueError as e:
            print(f"Error: Invalid RA/Dec format for custom target: {e}")
    elif json_data:
        # Iterate through constellations and objects to create FixedTarget instances
        for constellation, objects in json_data.items():
            for ra_str, dec_str, name in objects:
                try:
                    coords = ICRS(ra=ra_str, dec=dec_str) # Create ICRS coordinates
                    target = FixedTarget(coord=coords, name=name)
                    targets.append(target)
                except ValueError as e:
                    print(f"Warning: Could not parse coordinates for '{name}' ({ra_str}, {dec_str}): {e}. Skipping this target.")

    # Add the Moon as a dynamic target (its position changes with time)
    moon_coords = get_body('moon', times_utc, observer_location)
    moon_target = FixedTarget(coord=moon_coords, name="Moon")
    targets.append(moon_target)
    return targets

def plot_sky_darkness(ax, observer, times_utc_for_calc, times_plot_date):
    """
    Calculates and plots shaded regions representing different levels of sky darkness.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot on.
        observer (astroplan.Observer): The observer object.
        times_utc_for_calc (astropy.time.Time): A series of UTC times for sun altitude calculation.
        times_plot_date (astropy.time.Time): A series of times (e.g., local times) for plotting on x-axis.
    """
    # Define colors for different twilight/night phases
    color_civil = '#ADD8E6' # Light Blue
    color_nautical = '#6495ED' # Royal Blue
    color_astro_twilight = '#4169E1' # Dodger Blue
    color_astro_night = '#191970' # Midnight Blue

    # Calculate Sun's altitude for the given UTC times and observer location
    sun_altaz = observer.altaz(times_utc_for_calc, get_body('sun', times_utc_for_calc, observer.location))
    sun_alt = sun_altaz.alt.to(u.deg).value

    # Define altitude thresholds for twilight and night
    alt_astro_night = -18 # Sun 18 degrees below horizon (Astronomical Night)
    alt_astro_twilight_end = -12 # Sun 12 degrees below horizon (Astronomical Twilight ends)
    alt_nautical_twilight_end = -6 # Sun 6 degrees below horizon (Nautical Twilight ends)
    alt_civil_twilight_end = 0 # Sun at horizon (Civil Twilight ends)

    # Create masks for each darkness phase based on UTC calculated sun altitude
    mask_astro_night = sun_alt < alt_astro_night
    mask_astro_twilight = (sun_alt >= alt_astro_night) & (sun_alt < alt_astro_twilight_end)
    mask_nautical_twilight = (sun_alt >= alt_astro_twilight_end) & (sun_alt < alt_nautical_twilight_end)
    mask_civil_twilight = (sun_alt >= alt_nautical_twilight_end) & (sun_alt < alt_civil_twilight_end)

    # Helper function to plot shaded areas
    def plot_masked_areas(ax_target, plot_times, mask, color, label):
        # Fill between the current y-axis lower limit (or a bit below) and 90 degrees altitude
        ax_target.fill_between(plot_times.plot_date, ax_target.get_ylim()[0] - 30, 90,
                               where=mask, facecolor=color, alpha=0.5, label=label, zorder=0)

    # Plot the shaded areas for each darkness phase using the specified plot times
    plot_masked_areas(ax, times_plot_date, mask_astro_night, color_astro_night, 'Astronomical Night')
    plot_masked_areas(ax, times_plot_date, mask_astro_twilight, color_astro_twilight, 'Astronomical Twilight')
    plot_masked_areas(ax, times_plot_date, mask_nautical_twilight, color_nautical, 'Nautical Twilight')
    plot_masked_areas(ax, times_plot_date, mask_civil_twilight, color_civil, 'Civil Twilight')

def plot_object_altitudes(ax, observer, times_utc_for_calc, times_plot_date, targets):
    """
    Calculates and plots the altitude profiles for all celestial targets.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot on.
        observer (astroplan.Observer): The observer object.
        times_utc_for_calc (astropy.time.Time): A series of UTC times for altitude calculation.
        times_plot_date (astropy.time.Time): A series of times (e.g., local times) for plotting on x-axis.
        targets (list): A list of astroplan.FixedTarget objects.
    """
    for target in targets:
        # Calculate altitude and azimuth for each target at the given UTC times
        altaz_frames = observer.altaz(times_utc_for_calc, target)
        altitudes = altaz_frames.alt.to(u.deg).value

        # The target.name will already be the desired label (either from JSON or formatted RA/Dec)
        label_name = target.name

        # Plot the altitude data using the specified plot times for the x-axis
        if target.name == "Moon":
            ax.plot(times_plot_date.plot_date, altitudes, label=label_name, linestyle='--', zorder=1)
        else:
            ax.plot(times_plot_date.plot_date, altitudes, label=label_name, zorder=1)

def add_twilight_labels(ax, observer, plot_date, tomorrow_date, time_start_utc, time_end_utc, time_offset_local, times_local_for_plot):
    """
    Adds vertical lines and text labels to the plot indicating specific twilight start/end times.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot on.
        observer (astroplan.Observer): The observer object.
        plot_date (datetime.date): The date for the twilight calculations (evening).
        tomorrow_date (datetime.date): The date for the twilight calculations (morning).
        time_start_utc (astropy.time.Time): The start time of the plot in UTC.
        time_end_utc (astropy.time.Time): The end time of the plot in UTC.
        time_offset_local (astropy.units.quantity.Quantity): The local time offset from UTC.
        times_local_for_plot (astropy.time.Time): The local time array used for main x-axis plotting.

    Returns:
        dict: A dictionary containing the UTC Time objects for various twilight phases.
    """
    # Define colors for twilight labels
    color_civil = '#ADD8E6'
    color_nautical = '#6495ED'
    color_astro_twilight = '#4169E1'

    # Set up times for twilight calculations
    date_for_twilights_evening = Time(
        f"{plot_date.year}-"
        f"{plot_date.month:02d}-"
        f"{plot_date.day:02d}T00:00:00", # Start of the plot date for evening twilights
        format='isot', scale='utc'
    )
    date_for_twilights_morning = Time(
        f"{tomorrow_date.year}-"
        f"{tomorrow_date.month:02d}-{tomorrow_date.day:02d}T00:00:00", # Start of tomorrow for morning twilights
        format='isot', scale='utc'
    )

    # Calculate twilight start times (evening) in UTC
    twilight_start_civil_utc = observer.twilight_evening_civil(date_for_twilights_evening, which='next')
    twilight_start_nautical_utc = observer.twilight_evening_nautical(date_for_twilights_evening, which='next')
    twilight_start_astro_utc = observer.twilight_evening_astronomical(date_for_twilights_evening, which='next')

    # Calculate twilight end times (morning) in UTC
    twilight_end_astro_utc = observer.twilight_morning_astronomical(date_for_twilights_morning, which='nearest')
    twilight_end_nautical_utc = observer.twilight_morning_nautical(date_for_twilights_morning, which='nearest')
    twilight_end_civil_utc = observer.twilight_morning_civil(date_for_twilights_morning, which='nearest')

    def add_label(ax_target, time_utc, text_label, y_pos, color):
        """Helper function to add a vertical line and rotated text label."""
        # Only add label if the time (converted to local) falls within the plot's local time range
        # This ensures the vertical line is drawn in the correct local time position on the local axis
        if time_utc is not None: # Check if a twilight time was found
            time_local = time_utc + time_offset_local # Convert to local time for plotting and label
            if time_local.jd >= times_local_for_plot[0].jd and time_local.jd <= times_local_for_plot[-1].jd:
                ax_target.axvline(time_local.plot_date, color=color, linestyle=':', alpha=0.7, zorder=2)
                ax_target.text(time_local.plot_date, y_pos,
                               f"{time_local.strftime('%H:%M')} Local\n{text_label}",
                               rotation=90, va='bottom', ha='left', color=color, fontsize=8,
                               bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'), zorder=3)

    # Y-positions for the labels on the plot
    y_pos_civil = 3
    y_pos_nautical = 8
    y_pos_astro = 13

    # Add labels for evening twilights
    add_label(ax, twilight_start_civil_utc, "Civil Twilight Start", y_pos_civil, color_civil)
    add_label(ax, twilight_start_nautical_utc, "Nautical Twilight Start", y_pos_nautical, color_nautical)
    add_label(ax, twilight_start_astro_utc, "Astro Twilight Start", y_pos_astro, color_astro_twilight)

    # Add labels for morning twilights
    add_label(ax, twilight_end_astro_utc, "Astro Twilight End", y_pos_astro, color_astro_twilight)
    add_label(ax, twilight_end_nautical_utc, "Nautical Twilight End", y_pos_nautical, color_nautical)
    add_label(ax, twilight_end_civil_utc, "Civil Twilight End", y_pos_civil, color_civil)

    # Return the twilight times in UTC for use in the info box
    return {
        "civil_dusk": twilight_start_civil_utc,
        "civil_dawn": twilight_end_civil_utc,
        "nautical_dusk": twilight_start_nautical_utc,
        "nautical_dawn": twilight_end_nautical_utc,
        "astro_dusk": twilight_start_astro_utc,
        "astro_dawn": twilight_end_astro_utc,
    }


def add_twilight_info_box(ax, observer, twilight_times, time_offset_local):
    """
    Adds a text box to the top right corner of the plot summarizing twilight times.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot on.
        observer (astroplan.Observer): The observer object (added for sun_alt_min_utc calculation).
        twilight_times (dict): A dictionary containing UTC Time objects for various twilight phases.
        time_offset_local (astropy.units.quantity.Quantity): The local time offset from UTC.
    """
    twilight_info = "Twilight Times (Local):\n"

    def format_time_or_na(time_utc_obj):
        """Helper function to format a time object to HH:MM local, or 'N/A' if None."""
        if time_utc_obj is not None:
            return (time_utc_obj + time_offset_local).strftime('%H:%M')
        return "N/A"

    # Populate the info string with formatted twilight times
    twilight_info += f"  Civil Dusk: {format_time_or_na(twilight_times['civil_dusk'])}\n"
    twilight_info += f"  Civil Dawn: {format_time_or_na(twilight_times['civil_dawn'])}\n"
    twilight_info += f"  Nautical Dusk: {format_time_or_na(twilight_times['nautical_dusk'])}\n"
    twilight_info += f"  Nautical Dawn: {format_time_or_na(twilight_times['nautical_dawn'])}\n"

    # Only show astronomical twilight if it occurs (i.e., not a perpetual daytime scenario)
    astro_dusk_formatted = format_time_or_na(twilight_times['astro_dusk'])
    # This condition tries to check if astro twilight actually happens within a reasonable context
    # and if it's after civil dusk, preventing issues with 'nearest' or 'next' returning distant times.
    if twilight_times['civil_dusk'] and twilight_times['astro_dusk'] and (twilight_times['astro_dusk'].jd >= twilight_times['civil_dusk'].jd):
         twilight_info += f"  Astronomical Dusk: {astro_dusk_formatted}\n"
         twilight_info += f"  Astronomical Dawn: {format_time_or_na(twilight_times['astro_dawn'])}\n"

    # Check if any twilight occurred (more robust check, including potential "no night" scenarios)
    # This re-calculates sun altitude at civil dusk if available, otherwise defaults to 90 (sun always up)
    sun_alt_at_civil_dusk = observer.altaz(twilight_times['civil_dusk'], get_body('sun', twilight_times['civil_dusk'], observer.location)).alt.to(u.deg).value if twilight_times['civil_dusk'] else 90
    if sun_alt_at_civil_dusk > 0: # If sun never sets below 0 degrees at civil dusk
        twilight_info = "Twilight Times (Local):\n" # Reset info if no actual twilight
        twilight_info += "\n  No twilight or night occurs for this date/location\n  (Sun never sets below 0°)."

    # Add the text box to the plot
    ax.text(0.98, 0.98, twilight_info.strip(),
             transform=ax.transAxes, # Position relative to the axes
             fontsize=9,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='lightgray'))

def configure_plot(fig, ax1, ax2, plot_date, observer, time_offset_local, times_utc, times_local):
    """
    Configures the plot's titles, labels, legends, and axis formatting.

    Args:
        fig (matplotlib.figure.Figure): The main figure object.
        ax1 (matplotlib.axes.Axes): The primary matplotlib axes (bottom X-axis).
        ax2 (matplotlib.axes.Axes): The twin matplotlib axes (top X-axis).
        plot_date (datetime.date): The date used for the plot title.
        observer (astroplan.Observer): The observer object for location name.
        time_offset_local (astropy.units.quantity.Quantity): The local time offset from UTC.
        times_utc (astropy.time.Time): A series of UTC times for the top x-axis.
        times_local (astropy.time.Time): A series of local times for the bottom x-axis.
    """
    ax1.set_title(f"Altitude vs. Time at location {observer.name} ({plot_date.strftime('%Y-%m-%d')})")
    ax1.set_ylabel("Altitude (degrees)")
    ax1.grid(True)

    # Get handles and labels for the legend and sort them
    handles, labels = ax1.get_legend_handles_labels()
    unique_labels = {}
    for handle, label in zip(handles, labels):
        unique_labels[label] = handle

    # Custom sorting key to put darkness shading labels first
    def sort_key_darkness(label):
        if "Astronomical Night" in label: return 0
        if "Astronomical Twilight" in label: return 1
        if "Nautical Twilight" in label: return 2
        if "Civil Twilight" in label: return 3
        return 4 # For celestial objects, placed after twilight/night

    sorted_items = sorted(unique_labels.items(), key=lambda item: sort_key_darkness(item[0]))
    sorted_handles = [item[1] for item in sorted_items]
    sorted_labels = [item[0] for item in sorted_items]

    # Place the legend outside the plot area on the right
    ax1.legend(sorted_handles, sorted_labels, loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax1.set_ylim(-5, 90) # Set Y-axis limits from -5 to 90 degrees

    # Configure the bottom X-axis (ax1) for Local Time
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) # Only show hour:minute for local time
    ax1.set_xticks(times_local.plot_date[::len(times_local)//5]) # Use local times for ticks
    ax1.set_xlabel(f"Local Time ({observer.name}, UTC+{time_offset_local.to(u.hour).value:.0f}h)")


    # Configure the top X-axis (ax2) for UTC Time
    ax2.set_xlim(ax1.get_xlim()) # Ensure twin axis has same limits as primary
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) # Show hour:minute for UTC
    ax2.set_xticks(times_utc.plot_date[::len(times_utc)//5]) # Use UTC times for ticks
    ax2.set_xlabel(f"UTC Time ({plot_date.strftime('%Y-%m-%d')})")


    fig.autofmt_xdate() # Automatically format x-axis tick labels for better readability
    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for the legend

def main():
    """
    Main function to orchestrate the celestial observation plotting process.
    """
    plot_date, json_file_path, custom_ra, custom_dec = parse_arguments()
    # If date parsing failed or an invalid date was given that couldn't default to today, exit.
    if plot_date is None:
        return

    # Get observer, UTC times, Local times, and other time-related data
    observer, times_utc, times_local, time_start_utc, time_end_utc, time_offset_local, tomorrow_date = setup_observer_and_times(plot_date)

    if custom_ra and custom_dec:
        targets = create_targets(None, times_utc, observer.location, custom_ra, custom_dec)
    else:
        json_data = load_json_data(json_file_path)
        if json_data is None:
            return
        targets = create_targets(json_data, times_utc, observer.location)

    if not targets:
        print("No celestial objects to plot. Exiting.")
        return

    fig, ax1 = plt.subplots(figsize=(14, 8)) # Create figure and primary axes

    # Plot sky darkness and object altitudes using local times for the x-axis, but UTC for calculations
    plot_sky_darkness(ax1, observer, times_utc, times_local)
    plot_object_altitudes(ax1, observer, times_utc, times_local, targets)

    # Add twilight time labels and info box
    twilight_times = add_twilight_labels(ax1, observer, plot_date, tomorrow_date, time_start_utc, time_end_utc, time_offset_local, times_local)
    # Pass the observer object to add_twilight_info_box
    add_twilight_info_box(ax1, observer, twilight_times, time_offset_local)

    ax2 = ax1.twiny() # Create a twin x-axis for UTC time
    configure_plot(fig, ax1, ax2, plot_date, observer, time_offset_local, times_utc, times_local) # Configure plot aesthetics

    # Save and display the plot
    plt.savefig(f"celestial_altitudes_{plot_date.strftime('%Y-%m-%d')}.png", dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
