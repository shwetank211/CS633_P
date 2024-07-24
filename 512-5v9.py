import matplotlib.pyplot as plt
import numpy as np

# Data
data = {
    "512^2/": {
         "5": [0.059647,0.061515,0.061339,0.059412,0.051291]
    },
    "512^2//": {
        "9": [0.077007,0.077392,0.073359,0.074099,0.076102]
    }
}

# Extracting x and y data
x_data = []
y_data = []
stencil_colors = {'5': 'lightblue', '9': 'lightgreen'}
for size, values in data.items():
    for stencil, times in values.items():
        x_data.append(size + "-" + stencil)
        y_data.extend(times)

# Creating the boxplot
plt.figure(figsize=(12, 8))
bp = plt.boxplot([y_data[i:i+5] for i in range(0, len(y_data), 5)], positions=range(1, len(x_data) + 1), widths=0.6,
                 patch_artist=True)
plt.xticks(range(1, len(x_data) + 1), x_data, rotation=45)
plt.xlabel('Data Size and Stencil')
plt.ylabel('Time (in seconds)')
plt.title('Boxplot for Different Data Sizes and Stencils')
plt.grid(True)

# Customizing box colors
for box, color in zip(bp['boxes'], [stencil_colors[x.split('-')[1]] for x in x_data]):
    box.set(color='black', linewidth=1)
    box.set(facecolor=color)

# Add horizontal lines to connect median values
medians = [np.median(y_data[i:i+5]) for i in range(0, len(y_data), 5)]
for i, median in enumerate(medians):
    plt.plot([i + 0.5, i + 1.5], [median, median], color='gray', lw=2, ls='--')

# Add line graph connecting median values
for i in range(len(medians) - 1):
    plt.plot([i + 1.5, i + 2.5], [medians[i], medians[i + 1]], color='red', lw=2)


# Add labels to median values
for i, median in enumerate(medians):
    plt.text(i + 1, median, f'{median:.3f}', horizontalalignment='center', verticalalignment='bottom', fontsize=8)

plt.tight_layout()
plt.show()
