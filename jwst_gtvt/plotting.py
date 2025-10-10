from math import pi

from bokeh.io import output_file, save
from bokeh.layouts import gridplot
from bokeh.models import HoverTool, DatetimeTickFormatter
from bokeh.plotting import figure, show
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from astropy.time import Time

from jwst_gtvt.display_results import get_visibility_windows
from jwst_gtvt.utils import JWST_INSTRUMENTS

def plot_visibility(ephemeris, instrument=None, name=None, write_plot=None):
    """Make static visibility plot
    Parameters
    ----------
    ephemeris : jwst_gtvt.jwst_tvt.ephemeris
        Ephemeris class with fixed or moving target positions calculated.
    instrument : str
        JWST instrument name
    name : str
        Target name (designation from Horizons)
    write_plot : str
        Filename to write plot out to
    """
    dataframe = ephemeris.dataframe
    dataframe["times"] = Time(dataframe["MJD"], format="mjd").datetime

    df = dataframe.loc[dataframe["in_FOR"]]

    # These indices allow us to get the visible regions regions
    window_indices = get_visibility_windows(df.index.tolist())

    if instrument:
        # Set plotting configs
        plt.figure(figsize=(14, 8))
        plt.grid(color="k", linestyle="--", linewidth=2, alpha=0.3)
        plt.xticks(fontsize=14, rotation=45)
        plt.yticks(fontsize=14)

        for start, end in window_indices:
            data_to_plot = df.loc[start:end]
            min_PA_data = data_to_plot[instrument.upper() + "_min_pa_angle"]
            max_PA_data = data_to_plot[instrument.upper() + "_max_pa_angle"]
            plt.fill_between(
                data_to_plot["times"], min_PA_data, max_PA_data, color="grey"
            )
            plt.fmt_xdata = DateFormatter("%Y-%m-%d")

        if instrument.lower() == "v3pa":
            plt.ylabel(r"Available Position Angles ($^\circ$)", fontsize=18)
        else:
            plt.ylabel(r"Available Aperture Position Angles ($^\circ$)", fontsize=18)

        if ephemeris.fixed:
            ra, dec = max(df["ra"]), max(df["dec"])
            if name:
                plt.title("{} with {}".format(name, instrument.upper()), fontsize=18)
            else:
                plt.title(
                    "RA: {} Dec: {} with {}".format(
                        round(ra, 4), round(dec, 4), instrument.upper()
                    ),
                    fontsize=18,
                )
        else:
            plt.title(
                "Target {} with {}".format(ephemeris.target_name, instrument.upper()),
                fontsize=18,
            )

        if write_plot:
            plt.savefig(write_plot)
        else:
            plt.show()

    else:
        # plot all instruments here.
        fig, axs = plt.subplots(2, 3, figsize=(14, 8))

        if ephemeris.fixed:
            ra, dec = max(df["ra"]), max(df["dec"])
            if name:
                fig.suptitle("Target Name: {}".format(name), fontsize=18)
            else:
                fig.suptitle("RA: {} Dec: {}".format(ra, dec), fontsize=18)
        else:
            plt.suptitle("Target {}".format(ephemeris.target_name), fontsize=18)

        for instrument_name, ax in zip(JWST_INSTRUMENTS, axs.flatten()):
            for start, end in window_indices:
                data_to_plot = df.loc[start:end]
                min_PA_data = data_to_plot[instrument_name.upper() + "_min_pa_angle"]
                max_PA_data = data_to_plot[instrument_name.upper() + "_max_pa_angle"]
                ax.fill_between(
                    data_to_plot["times"], min_PA_data, max_PA_data, color="grey"
                )
                ax.fmt_xdata = DateFormatter("%Y-%m-%d")
                ax.set_title(instrument_name.upper())
                ax.tick_params("x", labelrotation=45)
                ax.grid(color="k", linestyle="--", linewidth=2, alpha=0.3)
                if instrument_name.lower() == "v3pa":
                    ax.set_ylabel("Available Position Angles (°)")
                else:
                    ax.set_ylabel("Available Aperture Position Angles (°)")

        fig.tight_layout()

        if write_plot:
            plt.savefig(write_plot)
        else:
            plt.show()


def plot_interactive_visibility(ephemeris, instrument=None, name=None, write_plot=None):
    """Make interactive visibility plot

    Parameters
    ----------
    ephemeris : jwst_gtvt.jwst_tvt.ephemeris
        Ephemeris class with fixed or moving target positions calculated.
    instrument : str
        JWST instrument name
    name : str
        Target name (designation from Horizons)
    write_plot : str
        Filename to write plot out to
    """

    def _make_plot(instrument, height, width, name=None):
        """Make bokeh plot with hover feature.

        Parameters
        ----------
        instrument : str
            JWST instrument name
        height : int
            Height size of plot in px
        width : int
            Width size of plot in px
        """
        if instrument.lower() == "v3pa":
            ylabel = "Available Position Angles (°)"
        else:
            ylabel = "Available Aperture Position Angles (°)"

        if name:
            title = f"{name} visibility for {instrument}"
        else:
            title = f"{instrument}"
        p = figure(
            title=title,
            x_axis_type="datetime",
            x_axis_label="Date",
            y_axis_label=ylabel,
            height=height,
            width=width,
        )
        for start, end in window_indices:
            data_to_plot = df.loc[start:end]

            p.varea(
                x="display_date",
                y1=f"{instrument.upper()}_min_pa_angle",
                y2=f"{instrument.upper()}_max_pa_angle",
                source=data_to_plot,
                fill_color="grey",
                fill_alpha=0.6,
            )

            # axis formatting
            p.xaxis.major_label_orientation = pi / 4
            p.title.text_font_style = "bold"
            p.title.text_font_size = "15pt"
            p.xaxis.axis_label_text_font_style = "bold"
            p.xaxis.axis_label_text_font_size = "15pt"
            p.yaxis.axis_label_text_font_style = "bold"
            p.yaxis.axis_label_text_font_size = "10pt"
            p.xaxis.major_label_text_font_style = "bold"
            p.xaxis.major_label_text_font_size = "12pt"
            p.yaxis.major_label_text_font_style = "bold"
            p.yaxis.major_label_text_font_size = "12pt"

            tooltips = [
                ("Date", "@display_date{%Y-%m-%d %H:%M:%S}"),
                ("Minimum PA", f"@{instrument.upper()}_min_pa_angle"),
                ("Maximum PA", f"@{instrument.upper()}_max_pa_angle"),
            ]
            formatters = {
                "@display_date": "datetime",
            }

            hover = HoverTool(
                tooltips=tooltips,
                formatters=formatters,
            )

            p.add_tools(hover)

        return p

    dataframe = ephemeris.dataframe

    df = dataframe.loc[dataframe["in_FOR"]]

    # These indices allow us to get the visible regions regions
    window_indices = get_visibility_windows(df.index.tolist())

    if instrument:
        single_instrument_plot = _make_plot(
            instrument, height=800, width=1200, name=name
        )
        if write_plot:
            output_file(write_plot)
            save(single_instrument_plot)
        else:
            show(single_instrument_plot)
    else:
        plots = []
        for instrument_name in JWST_INSTRUMENTS:
            plots.append(_make_plot(instrument_name, height=400, width=600, name=name))
        layout = gridplot(
            plots, ncols=3, merge_tools=False
        )  # Arrange in a grid with 3 columns
        if write_plot:
            output_file(write_plot)
            save(layout)
        else:
            show(layout)
