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
@INST= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
my $INST_index=3;
@DCA= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
my $DCA_index=4;
@DCR= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
my $DCR_index=5;
my $index = 0;
print "========= Threads=32 ==========\n";
print "Benchmark  DCRR  DCMR\n";
foreach (@benchmarks) {
  my $iterations=0;
  my $currbench = $_;
  print "$currbench";
  my $LOG = $LOGDIR . "/" . $currbench . ".3096.1024.pufferFish.threads-32.log";
  open(FILE, $LOG);
  my $flag = "false";
  foreach my $line (<FILE>) {
    chomp $line;
    if($line =~ /DATA_CACHE_ACCESSES/) {
      $flag = "true";
    } else {
      if($flag =~ /true/) {
        my @values = split('\t', $line);
	$INST[$index] += $values[$INST_index];
	$DCA[$index] += $values[$DCA_index];
	$DCR[$index] += $values[$DCR_index];
	$flag = "false";
	$iterations++;
      }
    }
  }
  $INST[$index] = $INST[$index] / $iterations;
  $DCA[$index] = $DCA[$index] / $iterations;
  $DCR[$index] = $DCR[$index] / $iterations;
  my $dcmr = $DCR[$index] / $INST[$index];
  my $dcrr = $DCA[$index] / $INST[$index];
  $dcmr = sprintf("%.4f",$dcmr);
  $dcrr = sprintf("%.2f",$dcrr);
  print "  $dcrr  $dcmr\n";
  close FILE;
  $index++;
}
