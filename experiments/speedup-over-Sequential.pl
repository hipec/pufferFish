#!/usr/bin/perl -w

#1) This script first calculates the average of "time.kernel" for
# the "sequential" execution of each benchmark.
#2) Then it calculate the average of "time.kernel" for each 
# benchmark, and for each of the builds, and for both threads 16 & 32.
#3) It then computes the speedup by dividing the sequential 
# time by parallel time for each benchmark and for each builds.

#If you changed the name of logs directory in experiment*.sh
#then change accordingly in both speedup calculation 
#scripts (including this)
my $LOGDIR = $ENV{PUFFERFISH} . "/logs";

@benchmarks=( "CilkSort", "SOR", "SORir", "LULESH", "LULESHir", "SRAD", "SRADir" );
@seqtime = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
@iterations = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

# Calculate sequential time first
my $index = 0;
foreach (@benchmarks) {
  my $currbench = $_;
  my $LOG = $LOGDIR . "/" . $currbench . ".3096.1024.sequential.threads-0.log";
  open(FILE, $LOG);
  my $flag = "false";
  $iterations[$index] = 0;
  foreach my $line (<FILE>) {
    chomp $line;
    if($line =~ /time.kernel/) {
      $flag = "true";
    } else {
      if($flag =~ /true/) {
        $iterations[$index]++;
        my @values = split('\t', $line);
	my $timek = $values[0];
	$seqtime[$index] += $timek;
	$flag = "false";
      }
    }
  }
  $seqtime[$index] = $seqtime[$index] / $iterations[$index];
  close FILE;
  $index++;
}

@threads=( 16, 32 );
@builds=( "defaultRand", "defaultRandiAll", "hptDA", "pufferFish", "cilk", "cilkiAll" );
@buildName=( "HClib(FT)", "HClib(IL)", "HClib(HPT_DA)", "PufferFish", "CilkPlus(FT)", "CilkPlus(IL)" );
foreach (@threads) {
my $currthread = $_;
print "================= Thread=$currthread =====================\n";
print "Benchmark";
foreach (@buildName) {
  print "  $_";
}
print "\n";
$index=0;
foreach (@benchmarks) {
  my $currbench = $_;
  my @results;
  print "$currbench  ";
  foreach (@builds) {
    my $currbuild = $_; 
    my $partime = 0.0;
    my $LOG=$LOGDIR."/".$currbench.".3096.1024.".$currbuild.".threads-".$currthread.".log";
    open(FILE, $LOG);
    my $flag = "false";
    my $total_iterations = 0;
    foreach my $line (<FILE>) {
      chomp $line;
      if($line =~ /time.kernel/) {
        $flag = "true";
      } else {
        if($flag =~ /true/) {
          $total_iterations++;
          my @values = split('\t', $line);
	  $partime += $values[0];
	  $flag = "false";
        }
      }
    }
    $partime = $partime / $total_iterations;
    close FILE;
    my $speedup = $seqtime[$index] / $partime;
    push(@results, $speedup);
  }
  $index++;
  foreach (@results) {
    my $rounded = sprintf("%.2f",$_);
    print "$rounded  ";
  }
  print "\n";
}
}

