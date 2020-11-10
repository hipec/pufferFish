#!/usr/bin/perl -w

#1) This script first calculates the average of "time.kernel" for
# the "sequential" execution of each benchmark.
#2) Then it calculate the average of "time.kernel" for each 
# benchmark, and for each of the builds, and for both threads 16 & 32.
#3) It then computes the ratio by dividing the sequential 
# time by parallel time for each benchmark and for each builds.

#If you changed the name of logs directory in experiment*.sh
#then change accordingly in both ratio calculation 
#scripts (including this)
my $LOGDIR = $ENV{PUFFERFISH} . "/logs-pcm";

@benchmarks=( "CilkSort", "SOR", "SORir", "LULESH", "LULESHir", "SRAD", "SRADir" );
#Col=0=> RAPL_PKG_ENERGY:PWR1
#Col=1=> L3_MISS:CPMC0
#Col=2=> L3_ACCESS:CPMC1
@pcmNames=( "Package Energy", "L3-cache Miss", "L2-cache Miss" );
@pcmActualNames=( "RAPL_PKG_ENERGY:PWR1", "L3_MISS:CPMC0", "L3_ACCESS:CPMC1" );

my $valueColumnIndex = 0;
foreach (@pcmActualNames) {
@hclibFT = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
# Calculate HClibFT first
my $index = 0;
foreach (@benchmarks) {
  my $iterations=0;
  my $currbench = $_;
  my $LOG = $LOGDIR . "/" . $currbench . ".3096.1024.defaultRand.threads-32.log";
  open(FILE, $LOG);
  my $flag = "false";
  foreach my $line (<FILE>) {
    chomp $line;
    if($line =~ /$pcmActualNames[$valueColumnIndex]/) {
      $flag = "true";
    } else {
      if($flag =~ /true/) {
        my @values = split('\t', $line);
	my $pcm = $values[$valueColumnIndex];
	$hclibFT[$index] += $pcm;
	$flag = "false";
	$iterations++;
      }
    }
  }
  $hclibFT[$index] = $hclibFT[$index] / $iterations;
  close FILE;
  $index++;
}

@builds=( "defaultRandiAll", "hptDA", "pufferFish" );
@buildName=( "HClib(IL)", "HClib(HPT_DA)", "PufferFish");
$index=0;
print "=========== $pcmNames[$valueColumnIndex] Relative to HCLIB (FT) ==========\n";
print "Benchmark";
foreach (@buildName) {
  print "  $_";
}
print "\n";
foreach (@benchmarks) {
  my $currbench = $_;
  my @results;
  print "$currbench  ";
  foreach (@builds) {
    my $currbuild = $_; 
    my $parvalue = 0.0;
    my $LOG=$LOGDIR."/".$currbench.".3096.1024.".$currbuild.".threads-32.log";
    open(FILE, $LOG);
    my $flag = "false";
    my $total_iterations = 0;
    foreach my $line (<FILE>) {
      chomp $line;
      if($line =~ /$pcmActualNames[$valueColumnIndex]/) {
        $flag = "true";
      } else {
        if($flag =~ /true/) {
          $total_iterations++;
          my @values = split('\t', $line);
	  $parvalue += $values[$valueColumnIndex];
	  $flag = "false";
        }
      }
    }
    $parvalue = $parvalue / $total_iterations;
    close FILE;
    my $ratio = $parvalue / $hclibFT[$index];
    push(@results, $ratio);
  }
  $index++;
  foreach (@results) {
    my $rounded = sprintf("%.2f",$_);
    print "$rounded  ";
  }
  print "\n";
}
$valueColumnIndex++;
}
