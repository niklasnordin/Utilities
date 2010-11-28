#!/usr/bin/perl

# Generate the cylinder blockMeshDict file for the von Karman LES case
#              16          17         18          19
#  y ^        0 - - - - - 1- - - - - 2 - - - - - 3
#    |        |           |          |           |  
#    --> x    |     0     |    1     |   2       |     ny1
#             |20         |21        |22         |23
#             4 - - - - - 5 - - - - -6 - - - - - 7
#             |           |          |           |  
#             |    3      |          |   4       |     ny2
#             |24         |25        |26         |27
#             8 - - - - - 9 - - - - 10 - - - - -11
#             |           |          |           |  
#             |    5      |   6      |   7       |     ny1
#             |28         |29        |30         |31
#             12 - - - - 13 - - - - 14 - - - - - 15
#                                                   
#                 nx1       nx2=ny2       nx3
use strict;

if ($#ARGV != 5)
{
    usage();
    exit;
}

my @input = @ARGV;
my $D = $ARGV[0];
my $H1 = $ARGV[1];
my $L1 = $ARGV[2];
my $L2 = $ARGV[3];
my $W = $ARGV[4];
my $dx = $ARGV[5];

my $w = 0.5*$W;
my $r = 0.5*$D;
my $r0 = 0.5*$D;
my $xe = $L2-$L1;

# aspect ratio for largest/smallest cells for the different blocks
my $g1 = 20.0;
my $g2 = 20.0;
my $g3 = 20.0;

my $nx1 = calcN($L1-$r, $g1, $dx);
my $nx2 = int($D/$dx) + 1;
my $nx3 = calcN($xe-$r, $g3, $dx);

my $ny1 = calcN($w-$r, $g2, $dx);
my $ny2 = $nx2;

my $nz = int(0.3*$H1/$dx) + 1;
#$nz = 1;

my $g1i = 1.0/$g1;
my $g2i = 1.0/$g2;
my $g3i = 1.0/$g3;

printHeader("2.0");
printConvertToMeters("1.0");
printVertices();
printBlocks();
printEdges();
printPatches();
printMerge();

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sub printHeader
{
    print "FoamFile\n";
    print "{\n";
    print "    version\t $_[0];\n";
    print "    format\t ascii;\n";
    print "    class\t dictionary;\n";
    print "    object\t blockMeshDict;\n";
    print "}\n\n";

    print "/*\n";
    print "Automatically generated with $0 @input\n";
    print "*/\n\n";
}

sub printConvertToMeters
{
    print "convertToMeters\t $_[0];\n\n";
}

sub printVertices
{
    print "vertices\n";
    print "(\n";
    print "    ( -$L1  $w 0 )\n";
    print "    ( -$r   $w 0 )\n";
    print "    (  $r   $w 0 )\n";
    print "    (  $xe  $w 0 )\n";

    print "    ( -$L1  $r 0 )\n";
    print "    ( -$r   $r 0 )\n";
    print "    (  $r   $r 0 )\n";
    print "    (  $xe  $r 0 )\n";

    print "    ( -$L1  -$r 0 )\n";
    print "    ( -$r   -$r 0 )\n";
    print "    (  $r   -$r 0 )\n";
    print "    (  $xe  -$r 0 )\n";

    print "    ( -$L1  -$w 0 )\n";
    print "    ( -$r   -$w 0 )\n";
    print "    (  $r   -$w 0 )\n";
    print "    (  $xe  -$w 0 )\n";

    print "    ( -$L1  $w $H1 )\n";
    print "    ( -$r   $w $H1 )\n";
    print "    (  $r   $w $H1 )\n";
    print "    (  $xe  $w $H1 )\n";

    print "    ( -$L1  $r $H1 )\n";
    print "    ( -$r   $r $H1 )\n";
    print "    (  $r   $r $H1 )\n";
    print "    (  $xe  $r $H1 )\n";

    print "    ( -$L1  -$r $H1 )\n";
    print "    ( -$r   -$r $H1 )\n";
    print "    (  $r   -$r $H1 )\n";
    print "    (  $xe  -$r $H1 )\n";

    print "    ( -$L1  -$w $H1 )\n";
    print "    ( -$r   -$w $H1 )\n";
    print "    (  $r   -$w $H1 )\n";
    print "    (  $xe  -$w $H1 )\n";

    print ");\n\n";
}

sub printBlocks
{
    print "blocks\n";
    print "(\n";
    print "    hex (4 5 1 0 20 21 17 16) ($nx1 $ny1 $nz) simpleGrading ($g1i $g2 1)\n";
    print "    hex (5 6 2 1 21 22 18 17) ($nx2 $ny1 $nz) simpleGrading (1 $g2 1)\n";
    print "    hex (6 7 3 2 22 23 19 18) ($nx3 $ny1 $nz) simpleGrading ($g3 $g2 1)\n";

    print "    hex (8 9 5 4 24 25 21 20) ($nx1 $ny2 $nz) simpleGrading ($g1i 1 1)\n";
    print "    hex (10 11 7 6 26 27 23 22) ($nx3 $ny2 $nz) simpleGrading ($g3 1 1)\n";

    print "    hex (12 13 9 8 28 29 25 24) ($nx1 $ny1 $nz) simpleGrading ($g1i $g2i 1)\n";
    print "    hex (13 14 10 9 29 30 26 25) ($nx2 $ny1 $nz) simpleGrading (1 $g2i 1)\n";
    print "    hex (14 15 11 10 30 31 27 26) ($nx3 $ny1 $nz) simpleGrading ($g3 $g2i 1)\n";

    print ");\n\n";
}

sub printEdges
{
    print "edges\n";
    print "(\n";
    print "/*\n";
    print "    arc 5 6 ( 0 $r0 0 )\n";
    print "    arc 6 10 ( $r0 0 0 )\n";
    print "    arc 9 10 ( 0 -$r0 0 )\n";
    print "    arc 9 5 ( -$r0 0 0 )\n";

    print "    arc 21 22 ( 0 $r0 $H1 )\n";
    print "    arc 26 22 ( $r0 0 $H1 )\n";
    print "    arc 25 26 ( 0 -$r0 $H1 )\n";
    print "    arc 25 21 ( -$r0 0 $H1 )\n";

    print "*/\n";  
    print ");\n\n";
}

sub printPatches
{
    print "patches\n";
    print "(\n";
    print "    wall walls\n";
    print "    (\n";

    print "        (5 21 22 6)\n";
    print "        (6 22 26 10)\n";
    print "        (10 26 25 9)\n";
    print "        (9 25 21 5)\n";

    print "    )\n\n";

    print "    patch inlet\n";
    print "    (\n";
    print "        (4 20 16 0)\n";
    print "        (8 24 20 4)\n";
    print "        (12 28 24 8)\n";
    print "    )\n";

    print "    patch outlet\n";
    print "    (\n";
    print "        (3 19 23 7)\n";
    print "        (7 23 27 11)\n";
    print "        (11 27 31 15)\n";
    print "    )\n";

    print "    wall sides\n";
    print "    (\n";
    print "        (16 17 1 0)\n";
    print "        (17 18 2 1)\n";
    print "        (18 19 3 2)\n";

    print "        (12 13 29 28)\n";
    print "        (13 14 30 29)\n";
    print "        (14 15 31 30)\n";

    print "    )\n";

    print "    wall lowerUpper\n";
    print "    (\n";

    print "        (0 1 5 4)\n";
    print "        (1 2 6 5)\n";
    print "        (2 3 7 6)\n";
    print "        (4 5 9 8)\n";
    print "        (6 7 11 10)\n";
    print "        (8 9 13 12)\n";
    print "        (9 10 14 13)\n";
    print "        (10 11 15 14)\n";

    print "        (20 21 17 16)\n";
    print "        (21 22 18 17)\n";
    print "        (22 23 19 18)\n";
    print "        (24 25 21 20)\n";
    print "        (26 27 23 22)\n";
    print "        (28 29 25 24)\n";
    print "        (29 30 26 25)\n";
    print "        (30 31 27 26)\n";

    print "    )\n";
    print ");\n\n";
}

sub printMerge
{
    print "mergPatchPairs\n";
    print "(\n";
    print ");\n";
}

sub calcN
{
    my $L = $_[0];
    my $R = $_[1];
    my $delta = $_[2];

    my $n = 1 + log($R)/log(($L-$delta)/($L-$R*$delta));
    return int($n);
}

sub usage
{
    print "usage:\n";
    print "\t $0 <D> <H1> <L1> <L2> <W> <delta>\n";
    print "\t where (recommended value):\n";
    print "\t D = diameter of cylinder. (0.04)\n";
    print "\t H1 = Height of the cylinder. (0.392)\n";
    print "\t L1 = Length from inlet to cylinder center. (0.4)\n";
    print "\t L2 = Length of computational domain. (1.36)\n";
    print "\t W = Width of the domain (0.56)\n";
    print "\t delta = cell size across the sides of the cylinder. (0.001)\n";
}
