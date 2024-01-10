import pandas as pd
import pylab
from sklearn.linear_model import LinearRegression
from tabulate import tabulate
from tkinter import *
from tkinter import filedialog as fd

# ################################# input #####################################

figure_rows = 1  # Output figure dimensions
figure_columns = 4
no_empty = True  # Don't plot empty wells
min_finalOD = 0.5
only_empty = False  # Only plot empty wells
max_finalOD = False
log_scale = True
dont_plot = ["Water", pylab.nan]  # List of sample names to ignore
title = "EK01 growth curves in various media 10^5 dilution"
target_od = 0.5
save = True # Save figure
calibrate = False

def calibration_curve(OD_measured):
    if OD_measured <= 0.33:
        OD_real = 3.38*OD_measured
    elif OD_measured > 0.33:
        OD_real = 3.38*0.33 + 6.8*(OD_measured - 0.33)
    else:
        OD_real = pylab.nan
    return OD_real

# ############################## Prepare data ##################################

root = Tk()
root.withdraw()
print("Select data.")
file = fd.askopenfilename()
root.destroy()
data = pd.read_excel(file)

data1 = data[list(data)[0]]
data2 = data[list(data)[1:]]

blank = data2.min(axis=1)

data2 = data2.sub(blank, axis=0)

data = pd.concat([data1, data2], axis=1)

z = list(dict.fromkeys([data.values[i][0] for i in range(len(data))]))
z = [i for i in z if i not in dont_plot]

if no_empty:
    data = data.drop(data2[data2.max(axis=1) < min_finalOD].index)

if only_empty:
    data = data.drop(data2[data2.max(axis=1) > max_finalOD].index)

data_copy = data.copy()

if calibrate:
    data1 = data[list(data)[0]]
    data2 = data[list(data)[1:]]
    data2 = data2.applymap(calibration_curve)
    data = pd.concat([data1, data2], axis=1)

mini = min([min(i) for i in data2.values if min(i)==min(i)])
maxi = max([max(i) for i in data2.values if max(i)==max(i)])

# ############################### Plot growth curve ####################################

x = list(data)[1:]

x_lab = list(data)[0]
y_lab = "ln(OD600)"

# Empty canvas in inches
fig = pylab.figure(figsize=(min(12, int(6*figure_columns/figure_rows)), min(6, int(12*figure_rows/figure_columns))+0.6))
# Shared X and Y axis labels
fig.supxlabel(x_lab, fontsize=18)
#fig.text(0.5, 0.02, x_lab, ha='center')
fig.supylabel(y_lab, fontsize=18)
#fig.text(0.01, 0.5, y_lab, va='center', rotation='vertical')
pylab.suptitle(title)
#fig.text(0.5, 0.96, title, ha='center', fontsize=13)

n = 1
# Per sample
for sample in z:
    # Dictates the destination of the next plot
    fig.add_subplot(figure_rows, figure_columns, n)
    n += 1

    # Per replicate
    replicates = data[data[list(data)[0]] == sample]
    for replicate in [replicates[i:i+1] for i in range(len(replicates))]:

        y = replicate.drop(list(data)[0], axis=1).values[0]

        if log_scale:
            pylab.scatter(x, [pylab.log(i) for i in y], s=3)
        else:
            pylab.scatter(x, y, s=3)
        # End of replicate

    x_lines = []
    for x_line in x_lines:
        pylab.axvline(x=x_line, c="black")
    pylab.title(sample)
    #pylab.ylim(mini-(maxi-mini)/40, maxi+(maxi-mini)/40)
    pylab.xlim(min(x), max(x)*1.05)
    # End of sample

# Prevent overlapping
pylab.tight_layout()
if save:
    pylab.savefig(title, dpi=300)
pylab.show()

# ############################## Calculate growth rate ############################

data = data_copy
y_lab = "Growth rate (1/h)"
title = "Growth rates"

# Empty canvas in inches
fig = pylab.figure(figsize=(min(12, int(6*figure_columns/figure_rows)), min(6, int(12*figure_rows/figure_columns))))
# Shared X and Y axis labels
fig.text(0.5, 0.02, x_lab, ha='center')
fig.text(0.01, 0.5, y_lab, va='center', rotation='vertical')
fig.text(0.5, 0.96, title, ha='center', fontsize=13)

# Final data init
results = []

n = 1
# Per sample
for sample in z:
    # Dictates the destination of the next plot
    fig.add_subplot(figure_rows, figure_columns, n)
    n += 1
    # Sample data init
    times = []
    rates = []

    # Per replicate
    replicates = data[data[list(data)[0]] == sample]
    for replicate in [replicates[i:i+1] for i in range(len(replicates))]:
        y = replicate.drop(list(data)[0], axis=1).values[0]
        # Replicate data init
        rate = []
        label = []
        time = 0

        # Per subset of 9 data points
        for i in range(len(y)-8):
            # Define subset
            window = y[i:i+9]
            # Filter subsets that contain OD measurements below 0.02
            if min(window) > 0.02:
                # Fit linear model
                reg = LinearRegression().fit(pylab.array(x[i:i+9]).reshape(-1, 1),
                                             pylab.array([pylab.log(j) for j in window]).reshape(-1, 1))
                c = reg.coef_[0][0]
                p = pylab.exp(reg.predict(pylab.array(x[i+4]).reshape(-1, 1)))
                label += "unlabeled",
            else:
                c = 0
                p = 0
                label += "low OD",
            if calibrate:
                p = calibration_curve(p)
            # Record time for target OD
            if time == 0 and p > target_od:
                time = x[i+4]
                times += time,
            # Record growth rate
            rate += c,
            # End of data

        rate_filtered = pylab.array(rate)[pylab.array(label) != "low OD"]
        exponential_limit = pylab.mean(sorted(rate_filtered, reverse=True)[:20])
        exponential_bool = (pylab.array(rate) > exponential_limit)*(pylab.array(label) != "low OD")
        index_min = min(exponential_bool.nonzero()[0])
        index_max = max(exponential_bool.nonzero()[0])
        while "low OD" in label[index_min:index_max+1]:
            index_min += 1 + min(exponential_bool[index_min+1:].nonzero()[0])
        label[index_min:index_max+1] = ["exponential"]*(index_max-index_min+1)
        label = ["post-exponential" if i == "unlabeled" else i for i in label]

        if max(y) >= min_finalOD:
            for lab in pylab.unique(label):
                i = pylab.where(pylab.array(label) == lab)[0]
                if lab == "exponential":
                    rates += pylab.mean(pylab.array(rate)[i]),
                    pylab.scatter(pylab.array(x[4:-4])[i], pylab.array(rate)[i], s=10, label=lab)
                else:
                    pylab.scatter(pylab.array(x[4:-4])[i], pylab.array(rate)[i], s=5, c="black")
        # End of replicate

    if len(replicates):
        results += (sample, round(float(pylab.mean(rates)), 2), round(float(pylab.median(times)), 1)),

    # Graph
    pylab.title(sample)
    pylab.ylim(0, 3.0)
    pylab.xlim(min(x), max(x))
    # End of sample

# Table
print(tabulate(results, headers=["Strain", "Growth rate (1/h)", f"Time to reach OD600 = {target_od} (h)"]))
# Graph
pylab.tight_layout(pad=2)
pylab.show()
