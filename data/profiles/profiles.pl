#!/usr/bin/perl

#
# profiles.pl  -  make generic profiles distributed with Scarse
#
# $Id: profiles.pl,v 1.1 2005/09/30 07:10:26 afrolov Exp $
#


################### Profile generator sequence #########################

my %pidx = ();		# Profile index
my $distr = "profiles";	# Profile distribution file

sub qjoin {		# Quote & join
	my $bra = shift; my $ket = shift; my $sep = shift;
	
	my @tmp = (); foreach (@_) { push @tmp, "$bra$_$ket"; }
	
	return join($sep, @tmp);
}

sub profile {		# Generate profile
	my ($file, $desc, $primaries, $illum, $gamma, $ipb) = @_;
	my $label = $desc; if (ref($desc) eq 'ARRAY') { ($desc, $label) = @$desc; }
	
	$ipb = "ipb -cd -M" unless $ipb;
	
	my $c = "$ipb";
	$c .= " -p$primaries" if $primaries;
	$c .= " -p$illum" if $illum;
	$c .= " -iRGB:$gamma" if $gamma;
	$c .= " -d'SCARSE: $desc'" if $desc;
	$c .= " -r'Copyright (C) 1999-2005 Scarse Project'";
	
	system "$c $file"; $pidx{$file} = $desc;
	
	return "<A HREF=\"$file\">$label</A>";
}

sub rgbspaces {		# Standard RGB color spaces
	my $dir = "spaces";
	
	mkdir $dir, 0755 if !(-d $dir);
	$pidx{"$dir/"} = "Standard RGB color spaces";
	
	return
	&qjoin("  <TR><TH>", "</TH></TR>", "\n",
		&profile("$dir/Adobe.icm", "Adobe RGB", 'Adobe', undef, undef),
		&profile("$dir/Apple.icm", "Apple RGB", 'Apple', undef, undef),
		&profile("$dir/ColorMatch.icm", "ColorMatch RGB", 'ColorMatch', undef, undef),
		&profile("$dir/ProPhoto.icm", "Kodak ProPhoto RGB", 'ProPhoto', undef, undef),
		&profile("$dir/sRGB.icm", "sRGB (simplified)", 'sRGB', undef, undef),
		&profile("$dir/WideGamut.icm", "Wide Gamut RGB", 'WideGamut', undef, undef)
	),
	&qjoin("  <TR><TH>", "</TH></TR>", "\n",
		&profile("$dir/Best.icm", "Best RGB", 'Best', undef, undef),
		&profile("$dir/Beta.icm", "Beta RGB", 'Beta', undef, undef),
		&profile("$dir/Bruce.icm", "Bruce RGB", 'Bruce', undef, undef),
		&profile("$dir/CIE.icm", "CIE RGB", 'CIE', undef, undef),
		&profile("$dir/Don4.icm", "Don RGB 4", 'Don4', undef, undef),
		&profile("$dir/ECI.icm", "ECI RGB", 'ECI', undef, undef),
		&profile("$dir/EktaSpace.icm", "Ekta Space RGB (PS5)", 'EktaSpace', undef, undef),
		&profile("$dir/NTSC.icm", "NTSC(1953) RGB", 'NTSC', undef, undef),
		&profile("$dir/PAL.icm", "PAL/SECAM RGB", 'PAL/SECAM', undef, undef),
		&profile("$dir/SMPTE-C.icm", "SMPTE-C(CCIR 601-1) RGB", 'SMPTE-C', undef, undef)
	);
}

sub mntrspaces {	# Generic monitor spaces
	my $dir = "display"; my @profiles = ();
	
	mkdir $dir, 0755 if !(-d $dir);
	$pidx{"$dir/"} = "Generic RGB monitor profiles";
	
	my %ps = (
		EBU		=> ['EBU/ITU', 'EBU'],
		HDTV		=> ['HDTV(CCIR 709)', 'HDTV'],
		P22		=> ['P22-EBU', 'P22'],
		Trinitron	=> ['Trinitron', 'TRIN']
	);
	
	foreach $p (sort keys %ps) {
		foreach $t (50, 55, 65, 75, 93) {
			my @mntr = ();
			
			foreach $g (1.5, 1.8, 2.2, 2.5) {
				my $name = "$ps{$p}[1]$t$g"; $name =~ s/\.//;
				my $desc = "Generic $ps{$p}[0] monitor; $t" . "00K, gamma $g";
				
				push @mntr, &profile("$dir/$name.icm", [$desc, $name], $p, "D$t", $g);
			}
			
			push @profiles, "<TH>$ps{$p}[0], $t" . "00K</TH>\n" . &qjoin("    <TD>", "</TD>", "\n", @mntr);
		}
	}
	
	return &qjoin("  <TR>", "\n  </TR>", "\n", @profiles);
}

sub makedist {		# Make distribution files
	my @dirs = sort(grep /\/$/, keys %pidx);
	
	system "tar czf $distr.tar.gz " . join(' ', @dirs);
	system "zip -q -r $distr.zip " . join(' ', @dirs);
}

sub makeidx {		# Make profile index
	open IDX, ">FILES";
	
	foreach (sort keys %pidx) { print IDX "$_\t $pidx{$_}\n"; }
	
	close IDX;
}



################### Generate profiles ##################################

my $date = `date +%m/%d/%Y`; chomp $date;
my $ipb =  `(ipb -h 2>&1) | head -1`; chomp $ipb;

my ($rgbprofiles, $addprofiles) = &rgbspaces();
my $mntrprofiles = &mntrspaces();

&makedist(); &makeidx();

my $tgz = `du -sk $distr.tar.gz`+0;
my $zip = `du -sk $distr.zip`+0;



################### HTML profile index #################################

print << "_END_HTML_";
<!-- \$Id\$ -->
<HTML><HEAD><TITLE>Scarse Profile Library</TITLE></HEAD>
<BODY BGCOLOR=white TEXT=black>

<H1>Scarse Project: Profile Library</H1>
(Updated $date; profiles generated with Scarse $ipb)

<HR>

<H3 ALIGN=center>Standard RGB color spaces:</H3>

<CENTER>
<TABLE BORDER=0 CELLPADDING=8 BGCOLOR=lavender>
$rgbprofiles
</TABLE>
</CENTER>

<H3 ALIGN=center>Additional RGB color spaces:</H3>

<CENTER>
<TABLE BORDER=0 CELLPADDING=8 BGCOLOR=lavender>
$addprofiles
</TABLE>
</CENTER>

<H3 ALIGN=center>Generic monitor profiles:</H3>

<CENTER>
<TABLE BORDER=0 CELLPADDING=8 BGCOLOR=lavender>
  <TR BGCOLOR=gray><TH>Phosphors, Color Temperature</TH>
    <TH>Gamma 1.5</TH><TH>Gamma 1.8</TH><TH>Gamma 2.2</TH><TH>Gamma 2.5</TH>
  </TR>
$mntrprofiles
</TABLE>
</CENTER>

<CENTER>
<P>Or download all profiles in
<A HREF="$distr.tar.gz">.tar.gz</A> ($tgz\kb) or
<A HREF="$distr.zip">.zip</A> ($zip\kb) format.
</CENTER>


<P ALIGN=center>[ Back to <A HREF="../../">Scarse Project Homepage</A> ]<HR>
<ADDRESS>Andrei Frolov &lt;<A HREF="mailto:andrei\@phys.ualberta.ca">andrei\@phys.ualberta.ca</A>&gt;</ADDRESS>

</BODY>
</HTML>
_END_HTML_

