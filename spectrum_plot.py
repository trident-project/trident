from matplotlib import pyplot

def plot_spectrum(wavelength, flux, filename):
    # number of rows and columns
    n_rows = 1
    n_columns = 1

    # blank space between edge of figure and active plot area
    top_buffer = 0.05
    bottom_buffer = 0.15
    left_buffer = 0.05
    right_buffer = 0.03
    
    # blank space between plots
    hor_buffer = 0.05
    vert_buffer = 0.05

    # calculate the height and width of each panel
    panel_width = ((1.0 - left_buffer - right_buffer - 
                    ((n_columns-1)*hor_buffer)) / n_columns)
    panel_height = ((1.0 - top_buffer - bottom_buffer - 
                     ((n_rows-1)*vert_buffer)) / n_rows)

    # create a figure (figsize is in inches)
    pyplot.figure(figsize=(12, 4))

    # loop over all panels
    for i in range(n_rows * n_columns):

        # get the row and column number
        my_row = i / n_columns
        my_column = i % n_columns

        # calculate the position of the bottom, left corner of this plot
        left_side = left_buffer + (my_column * panel_width) + \
          my_column * hor_buffer
        top_side = 1.0 - (top_buffer + (my_row * panel_height) + \
                          my_row * vert_buffer)
        bottom_side = top_side - panel_height

        # create an axes object on which we will make the plot
        my_axes = pyplot.axes((left_side, bottom_side, panel_width, panel_height))

        my_axes.plot(wavelength, flux)
        my_axes.xaxis.set_label_text("$\\lambda$ [$\\AA$]")
        my_axes.yaxis.set_label_text("relative flux")
        my_axes.set_xlim(wavelength.min(), wavelength.max())
        my_axes.set_ylim(0, 1.1)
        pyplot.savefig(filename)
