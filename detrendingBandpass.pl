#!/usr/bin/perl
use strict;
use List::Util qw( min max );

# AUTHOR:
# Scott Swanson
# Thomson Regenerative Biology Lab
# Morgridge Institute for Research
# April 04, 2019


mainSub();
exit;

# mainSub
# grabs/initializes the window sizes for the two filters
# calls buildLowPassWeights to generate gaussian filter weight matrix
# calls loadDataFromStdin
# iterates over samples, calling processSamples to perform bandpass filter
# iterates over time steps to print out rows of bandpass-filtered data
sub mainSub {
  my (@args)=(@ARGV);
  my $highPassWindow = $args[0] || 45;
  my $lowPassWindow = $args[1] || 3;

  my @lpWeight=();
  buildLowPassWeights($lowPassWindow, \@lpWeight);

  my $nTimeSteps=0;
  my @intensities=();
  my @tTimeStep=();
  $nTimeSteps = loadDataFromStdin(\@tTimeStep, \@intensities);

  my $nSamples = scalar @{$intensities[0]};
  my @filteredIntensities = ();
  for (my $iSample = 0; $iSample < $nSamples; $iSample++) {
    processSample(\@intensities, $iSample, $nTimeSteps, $highPassWindow, $lowPassWindow, \@lpWeight, \@filteredIntensities);
  }

  for (my $iTimeStep = 0; $iTimeStep < $nTimeSteps; $iTimeStep++) {
    print join("\t", $tTimeStep[$iTimeStep], @{$filteredIntensities[$iTimeStep]})."\n";
  }
  
}

# loadDataFromStdin
#   reads rows of intensity data from STDIN
#   stores results in two arrays whose addresses
#     have been passed as arguments
#   tTimeStep_a holds the clock times of the time steps
#   intensities_a_a is an array of arrays -- each row
#     is holds data for all samples for one time step
#     e.g. intensities_a_a->[$i][$j] is the intensity
#       for sample j at time step i.
# input data from STDIN
#   tab-delimited rows of numerical data
#   each row corresponds to a time step
#   each column corresponds to a sample
#     EXCEPT: First column contains clock time of the time step
#   any/all rows with '#' in the first character are echoed to
#     standard output without being parsed
#   the file is ASSUMED to be complete and correct.
sub loadDataFromStdin {
  my $tTimeStep_a = shift @_;
  my $intensities_a_a = shift @_;

  my $nTimeSteps = 0;
  while (<STDIN>) {
    if (! m/^#/) {
      chomp;
      my ($tTimeStep,@intensity)=split /\t/;
      $intensities_a_a->[$nTimeSteps]=[@intensity];
      $tTimeStep_a->[$nTimeSteps++]=$tTimeStep;
    }
    else {print}
  }

  return $nTimeSteps;
}

# processSample
# takes the input intensity arrays and runs the high pass and low pass
#   filters on all the time steps for one sample
# results are returned in the intensitiesAfterFilter_a_a array that
#   gets passed as an argument
sub processSample {
  my $intensities_a_a = shift @_;
  my $iSample = shift @_;
  my $nTimeSteps = shift @_;
  my $highPassWindow = shift @_;
  my $lowPassWindow = shift @_;
  my $lpWeight_a = shift @_;
  my $intensitiesAfterFilter_a_a = shift @_;

  my @intensitiesAfterHighPass = ();
  highPass($intensities_a_a, $iSample, $nTimeSteps, $highPassWindow, \@intensitiesAfterHighPass);
  lowPass(\@intensitiesAfterHighPass, $iSample, $nTimeSteps, $lowPassWindow, $lpWeight_a, $intensitiesAfterFilter_a_a);
}

# highPass
# calculates the result of running a high-pass filter over the input array
#   for one sample
# stores the results in the intensitiesAfterHighPass_a_a array passed as a reference
# the filter specification is:
#   for time step i
#     calculate the mean intensity for time steps [i-highPassWindow .. i+highPassWindow]
#     subtract the mean intensity from the original intensity at time step i
#   EDGE CASES
#     for time steps within highPassWindow steps of the beginning or end of the time course
#       calculate the mean only over the available region.
#     e.g., for a window of 20, for time step 7
#        calculate the mean intensity over time steps [0..27], and subtract from time step 7. 
sub highPass {
  my $intensities_a_a = shift @_;
  my $iSample = shift @_;
  my $nTimeSteps = shift @_;
  my $highPassWindow = shift @_;
  my $intensitiesAfterHighPass_a_a = shift @_;

  for (my $iTimeStep = 0; $iTimeStep < $nTimeSteps; $iTimeStep++) {
    my $iTimeStep0 = max(($iTimeStep-$highPassWindow), 0);
    my $iTimeStep1 = min( ($iTimeStep+$highPassWindow),($nTimeSteps-1));
    my $meanIntensity = 0;
    for my $iTimeStep ($iTimeStep0..$iTimeStep1) {
      $meanIntensity += $intensities_a_a->[$iTimeStep][$iSample];
    }
    $meanIntensity/=($iTimeStep1-$iTimeStep0+1);
    $intensitiesAfterHighPass_a_a->[$iTimeStep][$iSample]=$intensities_a_a->[$iTimeStep][$iSample]-($meanIntensity);
  }
}

# lowPass
# calculates the result of running a low-pass gaussian filter over the input array
#   for one sample
# gaussian weights are provided in lpWeight_a array.
# stores the results in the intensitiesAfterLowPass_a_a array passed as a reference
# the filter specification is:
#   for time step i
#     calculated the weighted sum of intensities over [i-lowPassWindow .. i+lowPassWindow]
#       using weights provided in the lpWeight_a array passed as an argument by reference.
#     calculate the mean intensity for time steps [i-highPassWindow .. i+highPassWindow]
#     subtract the mean intensity from the original intensity at time step i
#   EDGE CASES
#     for time steps within lowPassWindow steps of the beginning or end of the time course
#       calculate the mean only over the available region.
#     e.g., for a window of 3, for time step 1
#        calculate the weight sum intensity over time steps [0..4].
#     NOTE that this means values on the edges will be attenuated somewhat, since their
#        effective weighting matrix does not sum to 1.0. 
sub lowPass {
  my $intensities_a_a = shift @_;
  my $iSample = shift @_;
  my $nTimeSteps = shift @_;
  my $lowPassWindow = shift @_;
  my $lpWeight_a = shift @_;
  my $intensitiesAfterLowPass_a_a = shift @_;

  for (my $iTimeStep = 0; $iTimeStep < $nTimeSteps; $iTimeStep++) {
    my $nMinus= ($iTimeStep >= $lowPassWindow) ? $lowPassWindow : $iTimeStep;
    my $nPlus = (($iTimeStep+$lowPassWindow) < $nTimeSteps) ? $lowPassWindow : ($nTimeSteps-$iTimeStep-1);
    my $intensityWtdSum = 0;
    my $weightSum  = 0;
    if ($nMinus > 0) {
      for my $iTimeStepOffset (1..$nMinus) {
        $intensityWtdSum += $intensities_a_a->[$iTimeStep-$iTimeStepOffset][$iSample] * $lpWeight_a->[$iTimeStepOffset];
        $weightSum += $lpWeight_a->[$iTimeStepOffset];
      }
    }
    for my $iTimeStepOffset (0..$nPlus) {
      $intensityWtdSum += $intensities_a_a->[$iTimeStep+$iTimeStepOffset][$iSample] * $lpWeight_a->[$iTimeStepOffset];
      $weightSum += $lpWeight_a->[$iTimeStepOffset];
    }
    $intensitiesAfterLowPass_a_a->[$iTimeStep][$iSample]=sprintf("%.1f",($intensityWtdSum/$weightSum));
  }
}

# findVar
# fills the guassian weight matrix for the low pass filter
# first calls findVar to determine a good variance for
#   generating a weight matrix that sums to 1.0
# then calls normal to calculate the individual
#   values for each position in the weight matrix
# (yes, this duplicates work done by findVar, since it
#   calculates the corresponding values in order to
#   to find the sum of the weight matrix).
sub buildLowPassWeights {
  my $lowPassWindow = shift @_;
  my $lpWeight_a = shift @_;

  my $lpVariance = findVar($lowPassWindow);
  for my $iWeight (0..$lowPassWindow) {
     $lpWeight_a->[$iWeight]=normal(0,$lpVariance,$iWeight);
  }
}

# Iteratively determines a satisfactory variance
# for the low-pass filter's Gaussian weight matrix,
# given the width of the window, +/- n.
# (the objective is to choose a variance that produces
#  a sum over the weight matrix of very nearly 1.0)
# it begins by choosing an upper bound of n**2 and
#   a lower bound of 0.0
# it then begins a binary search between the lower
#   and upper bounds until the ratio of the two is
#   > 0.999999999
sub findVar {
  my $n=shift @_;

  my $done  = 0;
  my $upper  = $n**2;
  my $var = $upper / 2;
  my $lower = 0.0;
  my $sum = 0.0;
  my $preVar = $upper;
  while (! $done) {
    $sum = 0.0;
    for my $i (-$n..$n) {
      my $nor= normal(0,$var,$i);
      $sum+=$nor;
    }
    if ($sum >= 0.999) {
        $lower = $preVar = $var;
        $var = ($var + $upper)/2;
    }
    else {
        $upper = $preVar = $var;
        $var = ($var + $lower)/2;
    }
    if (0.999999999 < ($lower/$upper)){$done=1}
  }

  return $var;
}

# Normal Distribution values calculated as follows
# g(x)=( 1 / (sigma * sqrt( 2 * PI ))) * exp((-1/2) * ((x - mean)/sigma) ** 2)
# var == sigma**2
# PI/4 = atan2(1,1)
# 2*PI = 8 * atan2(1,1)
# g(x)=( 1 / sqrt(sigma**2 *  2 * PI )) * exp((-1/2) * ((x - mean) ** 2) / sigma**2)
# g(x)=( 1 / sqrt(var * 8 * atan2(1,1) )) * exp((-1/2) * ((x - mean) ** 2) / sigma**2)
# g(x)=( 1 / sqrt(8 * atan2(1,1) * var)) * exp((-1) * ((x - mean) ** 2) / (2 * var))
# g(x)= exp((-1) * ((x - mean) ** 2) / (2 * var)) / sqrt(8 * atan2(1,1) * var)
sub normal { 
  my ($mean, $var, $x) = (@_);
  return undef unless $var > 0; 
  return exp(-($x - $mean)**2/(2 * $var))  /   sqrt(8 * atan2(1,1) * $var);
} 
