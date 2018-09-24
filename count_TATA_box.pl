use Bio::SeqIO;
open (out,'>C:\Users\wbzheng\U.citrina\supply\contigs_v2_2telo_frame_modified_teloremoved.R200.count_TATA.10-15');
open (outb,'>C:\Users\wbzheng\U.citrina\supply\contigs_v2_2telo_frame_modified_teloremoved.R200.count_TATA.pos.10-15');
open (outc,'>C:\Users\wbzheng\U.citrina\supply\contigs_v2_2telo_frame_modified_teloremoved.R200.count_TATA.length_vs_number.10-15');
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => 'C:\Users\wbzheng\U.citrina\supply\contigs_v2_2telo_frame_modified_teloremoved.R200');
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$seq=$result_Trans->seq();
	$length=$result_Trans->length();
	@ds=split('',$seq);
	$former='A';
	$count=0;
	@cell=undef;
	$pos=0;
	foreach$keys(@ds){
		$this=$keys;
		$pos++;
		if(($this eq 'A' or $this eq 'T') and ($former eq 'A' or $former eq 'T')){
			$cell[$count]++;
			$ppp{$cell[$count]}=$length-$pos;#if it is 5 prime region, please use =$pos; 3 prime region, use =$length-$pos
			$former=$this;
		}
		if(($this eq 'A' or $this eq 'T') and ($former eq 'C' or $former eq 'G')){
			
			$count++;
			$cell[$count]++;
			$former=$this;
		}
		$former=$this;
	}
	@cell=sort{$b<=>$a}@cell;
	$nx{$cell[0]}++;
	$o3=$ppp{$cell[0]}+int($cell[0]/2); #if it is 5 prime region, please use -int***; 3 prime region, use +int***
	if($cell[0]>=10 and $cell[0]<15){#This is the limit of length
		$num_c{$o3}++;
	}
	print out "$idthis\t$cell[0]\t$o3";
	# foreach$keyb(@cell){
		# print out "$keyb\t";
	# }
	print out "\n";
}
foreach$keys(sort{$a<=>$b}keys%num_c){
	print outb "$keys\t$num_c{$keys}\n";
}
foreach$keys(sort{$a<=>$b}keys%nx){
	print outc "$keys\t$nx{$keys}\n";
}