#!/opt/local/bin/perl

#This version first anchors all A's then takes the highest scores and anchors B's that are in a range 0-75 away from A and then scores A by making
#sure the CPP comes after the KIP, thus eliminating the loophole of counting KIP and CPP as one.
use strict;
use Data::Dumper qw(Dumper);
use feature qw(say);
use List::Util qw( min max );
use Bio::Seq;
use Bio::SeqIO;

#motif A and B score thresholds, 30 max for motif A, 18 max for motif B.
my $inputA = 17;
my $inputB = 12;



print "File: ";
my $userinput = <STDIN>;
chomp $userinput;


#my $useroutput = "$userinput"."_out";
#chomp $useroutput;

my $seqio_obj = Bio::SeqIO->new(-file => "$userinput",
								-format => "fasta" 
								);
								
#my $seq_out = Bio::SeqIO->new(
#                              -file   => ">$useroutput",
#                             -format => "fasta",
#                             -width  => "32766"
#                              );
my %hashbest;
my %hashbestA;
my @blockA;
my @blockB;
my $elementsA;
my @allelementsA;
my $pB;
my $pA;
my $maxA;
my $elementsB;
my @allelementsB;
my @filteredA;
my @filteredB;
my $size;
my $sizeA;
my $GC;
my $proline;
my $maxA;
my $maxB;
my %hashA;
my @hashA;
my %hashB;
my %seen;
my @hashB;
my @positionA;
my @positionB;
my @difference;
my $score2c;
my $eachelementA;
my $eachkeyA;
my @eachelementA;
my @eachkeyA;
my $seqs;
my $keyA;
my $eachelementB;
my $eachkeyB;
my @eachelementB;
my @eachkeyB;
my $keyB;
my @seqs;     
my $posA;  
my $posB;
my $newB;
my @posA;
my @uniqB;
my @posB; 
my $trimA;
my $trimZ;
my $difference; 
my @posAtrim;  
my @posforA;
my @posforB;
my $score1;
my $scoreP;
my $scoreGC;
my $score2;
my $score3;
my $score4;
my $score5;
my $score6;
my $score7;
my $scorebacktrack;
my $score2b;
my $score1b;
my $score1a;
my $scoreB1;
my $scoreB2;
my $scoreB3;
my $scoreB4;
my $scoreB5;
my $scoreforC;
my $scoreforcaret;
my $scoreforend;
my $totalB;
my $total2;
my $uniqB;
my $allelementsB;
my @pB;
my $total1;
my $whatmatchedB;
my $whatmatchedA;
my $ids;
my @everything;
my $count12A;
my $length;
my $blue = "CTILMFWVA";
my @scoresforB;
my @scoresforA;
my $scorescan;
my $maxC;
my $maxbestb;
my $maxbestA;
my $locationofA;
my $locationofB;
my $scoreofA;
my $scoreofB;
                
while (my $seq_obj = $seqio_obj->next_seq)
{
	$seqs = $seq_obj->seq;	
	$ids = $seq_obj->id;	
	@seqs = ($seqs);	
	$length = $seq_obj->length;
			
				
		
	if ($length < 301)
	{	
		
				
				
		foreach $seqs (@seqs)	
		{	
			%hashA =();	
			%hashB =();	
			%seen =();
			$score3 =();	
			$count12A =0;	
			$newB = 0;
			$eachkeyA =0;
			$eachelementA =0;
			$scorescan = 0;
			@allelementsA =();	
			@allelementsB =();	
			@filteredA =();	
	     	@filteredB =();	
			@hashA =();	
			@hashB =();	
			@positionA =();	
			@positionB =();	
			@difference =();	
			@eachelementA =();	
			@eachkeyA =();	
			@eachelementB =();	
			@eachkeyB =();	
			@seqs =();     	
			@posA =();	
			@uniqB =();	
			@posB =(); 	
			@posforA =();
			@posforB =();
			@pB =();	
			@everything =();	
			@scoresforB =();	
			@scoresforA =();
			@blockA =();
			@blockB =();	
			%hashbest =();
			%hashbestA =();
			$maxbestb = 0;
			$maxbestA = 0;
			$locationofA = 0;
			$locationofB = 0;
			$scoreofA = 0;
			$scoreofB = 0;
			
			while ($seqs =~ /[DEG][DECG]/g)				{push @positionA, $-[0];} 	
			while ($seqs =~ /[YCILMFWVA][CILMFWVA]/g)	{push @positionB, $-[0] - 12;}	
			
			sub uniq
			{%seen;	grep !$seen{$_}++, @_;}	
#--------------------------------------------------------------------------------------------------------------#	
#					         		Grab Substring for A													   #
#--------------------------------------------------------------------------------------------------------------#	
					
			foreach (@positionA)	
			{	
				push @allelementsA, (substr $seqs, $_, 30);	
			}	
								
			@hashA{@positionA} = @allelementsA;	
			#print Dumper \%hashA;	
			
			
#--------------------------------------------------------------------------------------------------------------#	
#											MOTIF A															   #
#--------------------------------------------------------------------------------------------------------------#	
			
			
			$eachkeyA = keys %hashA;	
		   	$eachelementA = values %hashA;	
					
			while(($eachkeyA, $eachelementA) = each(%hashA))
			{	
			   	$score1 = 0; $score2 = 0; $score3 = 0; $score4 = 0; $score5 =0; $score1a = 0; $score1b = 0; $score2b = 0; $score2c = 0; $scoreforC = 0; $scoreforcaret = 0; $scoreforend = 0; 
			
				if ($eachelementA =~ /^[DE][DE].*C.?TP/)
				{$scoreforcaret = 4;
				 $scoreforC = 2;
						if 	  ($eachelementA =~ /TP.{0,8}?[KR][$blue]P/){$score1 = 10;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
											
						elsif ($eachelementA =~ /TP.{0,8}?[KR][$blue]/){$score2 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
									
						elsif ($eachelementA =~ /TP.{0,8}?[KR].P/){$score3 = 8;
								if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
								
						elsif ($eachelementA =~ /TP.{0,8}?[$blue]P/){$score4 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						elsif ($eachelementA =~ /TP.{0,8}?[RK]/){$score5 = 6;	
								if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0; $score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						
						
						
						else{my $cumscore1 = 0;}
				
						
						
					
			   	}	
			   	
			   	

                #------------------------------------------------------------------------------------#	
			   		
				elsif ($eachelementA =~ /^[DE][DE].*TP/){
				$scoreforcaret = 4;
				$scoreforC = 0;
						if 	  ($eachelementA =~ /TP.{0,8}?[KR][$blue]P/){$score1 = 10;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
											
						elsif ($eachelementA =~ /TP.{0,8}?[KR][$blue]/){$score2 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
									
						elsif ($eachelementA =~ /TP.{0,8}?[KR].P/){$score3 = 8;
								if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
								
						elsif ($eachelementA =~ /TP.{0,8}?[$blue]P/){$score4 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						elsif ($eachelementA =~ /TP.{0,8}?[RK]/){$score5 = 6;	
								if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0; $score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						
						
						
						else{my $cumscore1 = 0;}
				
						
						
					
				}	
				
				
                #------------------------------------------------------------------------------------#	
			   		
				elsif ($eachelementA =~ /C.?TP/){
				$scoreforcaret = 0;
				$scoreforC = 2;
						if 	  ($eachelementA =~ /TP.{0,8}?[KR][$blue]P/){$score1 = 10;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
											
						elsif ($eachelementA =~ /TP.{0,8}?[KR][$blue]/){$score2 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
									
						elsif ($eachelementA =~ /TP.{0,8}?[KR].P/){$score3 = 8;
								if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
								
						elsif ($eachelementA =~ /TP.{0,8}?[$blue]P/){$score4 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						elsif ($eachelementA =~ /TP.{0,8}?[RK]/){$score5 = 6;	
								if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0; $score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						
						
						
						else{my $cumscore1 = 0;}
				
						
						
						
				}	
				
                #------------------------------------------------------------------------------------#	
			   		
				elsif ($eachelementA =~ /TP/){
				$scoreforcaret = 0;
				$scoreforC = 0;
						if 	  ($eachelementA =~ /TP.{0,8}?[KR][$blue]P/){$score1 = 10;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
											
						elsif ($eachelementA =~ /TP.{0,8}?[KR][$blue]/){$score2 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;
									if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR][$blue].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
									
						elsif ($eachelementA =~ /TP.{0,8}?[KR].P/){$score3 = 8;
								if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[KR].P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
								
						elsif ($eachelementA =~ /TP.{0,8}?[$blue]P/){$score4 = 8;	
								if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0;	$score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[$blue]P.{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						elsif ($eachelementA =~ /TP.{0,8}?[RK]/){$score5 = 6;	
								if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP/){$score1a = 5;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P/){$score1b = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR]/){$score2b = 2;
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR]/){$score2c = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}CP.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
												else{$scoreforend = 0;}}
											else{$score2c = 0;}}
										else{$score2b = 0;}
									}	
									else{$score1b = 0;}
								}		
								else{$score1a = 0; $score1b = 0; $score2b = 0;	
									if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P/){$score1a = 3;	
										if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P/){$score1b = 3;	
											if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR]/){$score2b = 2;
												if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR]/){$score2c = 2;
													if ($eachelementA =~ /TP.{0,8}?[RK].{1,20}[PST]P.?.P.?[KR][KR].{0,2}[PRKQ]/){$scoreforend = 2;}
													else{$scoreforend = 0;}}
												else{$score2c = 0;}}
											else{$score2b = 0;}
										}
										else{$score1b = 0;}
									}
									else{$score1a = 0;}
								}
						}
						
						
						
						
						else{my $cumscore1 = 0;}
				
						
						
					
				}	
				
							
				
				
				
							
				$total1 = $score1 + $score2 + $score3 + $score4 + $score5;	
				$total2 = $score1a + $score1b + $score2b + $score2c + $scoreforC + $scoreforcaret + $scoreforend;
				my $cumscore1 = ($total1 + $total2);
				
					
					if (1){	
					#print ">$ids|$total1|$total2|=$cumscore1\n$eachelementA\n";	
					push @posforA, $eachkeyA;	
					push @scoresforA, $cumscore1;
					push @blockA, $eachelementA;	
					}	
				}	#say Dumper \@posforA;
					#say Dumper \@scoresforA;
				


#--------------------------------------------------------------------------------------------------------------#	
#											MOTIF B															   #
#--------------------------------------------------------------------------------------------------------------#	
		foreach $pB (@positionB)	{	
				foreach $pA (@posforA) {	
					if ($pB > $pA+30){push @pB, $pB;}}}	
						
					@pB = uniq(@pB);
					#say Dumper \@pB;

					foreach (@pB){	
					push @allelementsB, (substr $seqs, $_, 29);}	
			
					my %hashB;	
					@hashB{@pB} = @allelementsB;	
			
					$eachkeyB = keys %hashB;	
		   			my $eachelementB = values %hashB;	
		   			#say Dumper \%hashB;
		   				
					while(($pB, $eachelementB) = each(%hashB)) {	
						$scoreB1 = 0; $scoreB2 = 0; $scoreB3 = 0; $scoreB4 = 0; $scoreB5 = 0; $scorebacktrack = 0;	
							
						if ($eachelementB =~ /[RKQ][RKQ][RKQ].{0,15}[DE][CILMFWVADE][DE]/){$scorebacktrack = 6;
							if 	  ($eachelementB =~ /[DE][CILMFWVADE][DE]/){$scoreB1 = 8;	
								if ($eachelementB =~ /[DE][CILMFWVADE][DE].?[CILMFWVA]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE][CILMFWVADE][DE].?[CILMFWVA][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}	
							elsif ($eachelementB =~ /[DE][CILMFWVADE]/){$scoreB2 = 5;	
								if ($eachelementB =~ /[DE][CILMFWVADE].?[CILMFWVA]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE][CILMFWVADE].?[CILMFWVA][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}	
							elsif ($eachelementB =~ /[DE]/){$scoreB3 = 3;	
								if ($eachelementB =~ /[DE].?[CILMFWVADE]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE].?[CILMFWVADE][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}}
								
								
								
						elsif ($eachelementB =~ /[RKQ][RKQ].{0,15}[DE][CILMFWVADE][DE]/){$scorebacktrack = 4;
							if 	  ($eachelementB =~ /[DE][CILMFWVADE][DE]/){$scoreB1 = 8;	
								if ($eachelementB =~ /[DE][CILMFWVADE][DE].?[CILMFWVA]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE][CILMFWVADE][DE].?[CILMFWVA][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}	
							elsif ($eachelementB =~ /[DE][CILMFWVADE]/){$scoreB2 = 5;	
								if ($eachelementB =~ /[DE][CILMFWVADE].?[CILMFWVA]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE][CILMFWVADE].?[CILMFWVA][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}	
							elsif ($eachelementB =~ /[DE]/){$scoreB3 = 3;	
								if ($eachelementB =~ /[DE].?[CILMFWVADE]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE].?[CILMFWVADE][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}}
								
								
						
						else  {$scoreB1 = 0; $scoreB2 = 0; $scoreB3 = 0; $scoreB4 = 0; $scoreB5 = 0; $scorebacktrack = 0;
							if 	  ($eachelementB =~ /[DE][CILMFWVADE][DE]/){$scoreB1 = 8;	
								if ($eachelementB =~ /[DE][CILMFWVADE][DE].?[CILMFWVA]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE][CILMFWVADE][DE].?[CILMFWVA][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}	
							elsif ($eachelementB =~ /[DE][CILMFWVADE]/){$scoreB2 = 5;	
								if ($eachelementB =~ /[DE][CILMFWVADE].?[CILMFWVA]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE][CILMFWVADE].?[CILMFWVA][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}	
							elsif ($eachelementB =~ /[DE]/){$scoreB3 = 3;	
								if ($eachelementB =~ /[DE].?[CILMFWVADE]/){$scoreB4 = 2;}	
								else{$scoreB4 = 0;}	
								if ($eachelementB =~ /[DE].?[CILMFWVADE][CILMFWVA]/){$scoreB5 = 2;}	
								else{$scoreB5 = 0;}}}	
						
						$size = length($eachelementB);	
						
				my $totalB = ($scoreB1 + $scoreB2 + $scoreB3 + $scoreB4 + $scoreB5 + $scorebacktrack);	
				#say $totalB;
				if (1){
					push @posforB, $pB;
					push @scoresforB, $totalB;
					push @blockB, $eachelementB;
					}


				}

			}				

						
						my $maxA = max @scoresforA;	
						my $maxB = max @scoresforB;	
						my $AandB = $maxA + $maxB;
						#say $maxB;
						
						##################MOTIF A SCORING########################
						my %hashbestA;
						@hashbestA{@posforA} = @scoresforA;			
						my $maxbestA = max(values %hashbestA);
						my %hash_maxA = map { $hashbestA{$_}==$maxbestA ? ($_, $maxbestA) : () } keys %hashbestA;
						
						#print Dumper \%hashbestA;
						#print Dumper \%hash_maxA;

						while (my ($keyA, $valueA) = each (%hash_maxA))
		   				{
		   				#say ">$ids"."|Motif_A|Score: "."$valueA"."/30";
		   				#last;
						#$seqs = lc $seqs;
						#substr( $seqs, $_, 1 ) =~ tr[a-z][A-Z] for $keyA..$keyA+30;
						#print "$seqs\n\n";
		   				}
						########################################################
						
						#say Dumper \@posforB;
						#say Dumper \@scoresforB;
						###################Motif B Scoring########################
						my %hashbest;	
						@hashbest{@posforB} = @scoresforB;
						#say Dumper \%hashbest;
						my $maxbestb = max(values %hashbest);
						my %hash_max = map { $hashbest{$_}==$maxbestb ? ($_, $maxbestb) : () } keys %hashbest;
	
						#print Dumper \%hashbest;
						#print Dumper \%hash_max;

						while (my ($key, $value) = each (%hash_max))
		   				{
		   				#say ">$ids"."|Motif_B|Score: "."$value"."/18";
						#$seqs = lc $seqs;
						#substr( $seqs, $_, 1 ) =~ tr[a-z][A-Z] for $key..$key+29;;
						#say $seqs;;
		   				#last;
		   				}
						#########################################################
						
						
						
						#my ($keyA, $valueA) = each (%hash_maxA))
						
						while (($locationofA, $scoreofA) = each (%hash_maxA))
						{
							while (($locationofB, $scoreofB) = each (%hash_max))
							{
								if (($scoreofA >= $inputA) && ($scoreofB >= $inputB) && ($locationofB > $locationofA+30))
								{
								my $sum = $scoreofA+$scoreofB;
								say ">$ids"."|Motif_A/B_Score:"."$scoreofA"."/30"."_&_"."$scoreofB"."/18"."|$sum|";
								$seqs = lc $seqs;
								substr( $seqs, $_, 1 ) =~ tr[a-z][A-Z] for $locationofA..$locationofA+30,$locationofB..$locationofB+29;
								say $seqs;
								}
								
							}

						}

						
	}

}
