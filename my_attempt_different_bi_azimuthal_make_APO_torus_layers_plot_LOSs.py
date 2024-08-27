import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import matplotlib as mpl

jupiter_image = Image.open('jpegPIA07783.png')
#jupiter_image = jupiter_image.crop((200, 200, 1400, 1400))  # Crop to remove black edges
#jupiter_image = jupiter_image.resize((200, 200))  # Resize to fit the plot



# Define the radii
#radii = [4.5, 5., 5.5, 6., 6.5, 7., 7.5]
radii = [ 4.5, 5., 5.5, 6., 6.5, 7., 7.5,8.]
#additional_radius = 8
vertical_lines_x = [7.75, 7.25, 6.75, 6.25, 5.75, 5.25, 4.75, -7.75, -7.25, -6.75, -6.25, -5.75, -5.25, -4.75]

# Define colors using the updated method
#colors_left = plt.colormaps['tab10'](np.linspace(0, 1, 7))
#colors_right = plt.colormaps['viridis'](np.linspace(0, 1, 7))

# Define colors using the updated method
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=1, vmax=14)

colors = cmap(norm(np.linspace(1, 14, 14)))
colors_left = colors[0:7, :]
colors_right =colors[7:, :]

# Create a figure and axis
fig, ax = plt.subplots()

# Add grid for better visualization
ax.grid(True)


# Create 1D arrays for x and y
x = np.arange(-8, 8.01, 0.01)
y = np.arange(-8, 8.01, 0.01)

# Create a meshgrid for x and y
X, Y = np.meshgrid(x, y)

# Calculate the radius for each point
R = np.sqrt(X**2 + Y**2)

# Initialize the z array with NaNs
Z = np.full(X.shape, np.nan)
val=4
# Define the conditions and values based on the problem statement
conditions = [
    (X < -val) & (R > 4.5) & (R < 5.),
    (X < -val) & (R > 5.0) & (R < 5.5),
    (X < -val) & (R > 5.5) & (R < 6.0),
    (X < -val) & (R > 6.0) & (R < 6.5),
    (X < -val) & (R > 6.5) & (R < 7.0),
    (X < -val) & (R > 7.0) & (R < 7.5),
    (X < -val) & (R > 7.5) & (R < 8.0),
    (X >= val) & (R > 4.5) & (R < 5.),
    (X >= val) & (R > 5.0) & (R < 5.5),
    (X >= val) & (R > 5.5) & (R < 6.0),
    (X >= val) & (R > 6.0) & (R < 6.5),
    (X >= val) & (R > 6.5) & (R < 7.0),
    (X >= val) & (R > 7.0) & (R < 7.5),
    (X >= val) & (R > 7.5) & (R < 8.0)
]

# Values for each condition
values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

# Apply conditions and values to the Z array
for cond, value in zip(conditions, values):
    Z[cond] = value

ax.contourf(X, Y, Z, levels=15, cmap='viridis')



# Plot each circle
for radius in radii:
    circle = plt.Circle((0, 0), radius, fill=False)
    ax.add_artist(circle)
# Plot each vertical line
for x_line in vertical_lines_x:
    ax.axvline(x_line, color='grey', linestyle='--')


# Add labels for Dawn and Dusk at specified positions
ax.text(-6.25, 8.25, 'Dawn', fontsize=12, verticalalignment='top', horizontalalignment='center', color='black')
ax.text(6.25, 8.25, 'Dusk', fontsize=12, verticalalignment='top', horizontalalignment='center', color='black')

# Set the limits of the plot
ax.set_xlim(-9, 9)
ax.set_ylim(-9, 9)

# Customize tick locations and labels
ax.set_xticks(np.arange(-8, 9, 2))
ax.set_yticks(np.arange(-8, 9, 2))

# Ensure the aspect ratio is equal
ax.set_aspect('equal')



# Add the image of Jupiter's north pole at the center
image_extent = [-1, 1, -1, 1]
ax.imshow(jupiter_image, extent=image_extent, aspect='equal', zorder=5)


# Show the plot with updated labels and title
plt.title("Io Plasma Torus Lines of Sight")
plt.xlabel(r"x ($R_J$)")
plt.ylabel(r"y ($R_J$)")
#plt.set_aspect('equal')
plt.savefig('my_attempt_different_APO_LOSs_visualizations_final.png', dpi=1000)
plt.show()
