
# Detrending Bandpass 
detrendingBandpass.pl implements a very simple bandpass filter on a time-series dataset

It is intended for "detrending" and smoothing data to help visualize and even quantify periodic signals, typically in datasets where FFT is not useful. (For example, when the frequency is not constant.) It was originally developed for analysis of oscillating signals in luminescence data taken from cell culture plates. In this context, the oscillations can be hidden by the slow degradation/consumption of the luminescing reagents, and the signal can be relatively noisy due to the gross scale of the optical measurements.

The script contains no error-checking, nor complex arguments-handling.
The script takes input only from STDIN.
The script writes output only to STDOUT (and STDERR).

## Installing / Getting started

The script should run as written, in any environment with a Perl 5 installation. No installation is necessary.
To test the script, run it on the testInput.tab file in the testData directory, and use "diff" to compare the results with the provided expectedOutput.tab file.
Diff should report no differences.

```shell
git clone https://github.com/stewart-lab/Detrending.git
cd testData
../detrendingBandpass.pl < testInput.tab > filteredOutput.tab
diff filteredOutput.tab expectedOutput.tab | head
```

... should output nothing at all.

## Features

### Input
The algorithm processes one or more time series, input as a tab-delimited text matrix of values, where:
* Each row is a series.
* The first column contains a label for that series.
* All series have the same number of time steps.
* All time steps within a series are equal in duration.

### Output
* Any input row beginning with '#' is echoed directly to the output.
* Each filtered time series is output as a tab-delimited row, corresponding 1 to 1 with the original input file, including the initial series-label column

## Configuration

The script accepts only two numerical arguments from the command-line: The high-pass window and the low-pass window (See below for full explanations.) These must be provided as integers, in order. If only one is present, it will be processed as the high-pass window.

#### highPassWindow
Type: `Decimal integer`  
Default: `45`

For each time point in a series, the filter calculates the mean over a window of +/- highPassWindow neighboring time points, and subtracts that mean from the original value of the time point.

For time points near the edges of the series (ie. < highPassWindow from the edge), the window is accordingly truncated on that side.

Example:
```bash
./detrendingBandpass.pl 30 < myTimeSeries.tab > myFilteredTimeSeries_30_3.tab
```

Filters each time series in myTimeSeries.tab using a high-pass window of +/- 30, and a low-pass window of +/- 3 (the default).

#### Low-pass window-size
Type: `Decimal integer`  
Default: `3`

Can only be used if the high-pass window-size is also explicitly set on the command line.
For each time point in a series, the filter calculates a weighted Gaussian sum over a window of +/- lowPassWindow neighboring time points.

For time points near the edges of the series (ie. < lowPassWindow from the edge), the window is accordingly truncated on that side. THIS WILL ATTENUATE
THE SIGNAL NEAR THE EDGES OF THE SERIES. The gaussian weight-matrix is calculated dynamically during initialization. The weights are chosen to sum to 1.0 across the entire window. No adjustment is made for time steps at the periphery, so the sum of the weights in those cases will be less than 1.0, resulting in an overall attenuation of the the signal.

Example:
```bash
./detrendingBandpass.pl 25 5 < myTimeSeries.tab > myFilteredTimeSeries_25_5.tab
```

Filters each time series in myTimeSeries.tab using a high-pass window of +/- 25, and a low-pass window of +/- 5.

## Links

- Project homepage: https://github.com/stewart-lab/Detrending/
- Repository: https://github.com/stewart-lab/Detrending/
- Issue tracker: https://github.com/stewart-lab/Detrending/issues

## Licensing

