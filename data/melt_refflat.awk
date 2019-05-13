# Melts a refFlat file by decomposing exonstarts, exonends, and exonframes.
# An example of refflat can be obtained via:
#  

BEGIN { OFS="\t" }

NR==1 {
  print
}

NR>1 {
  # construct columns 1-9 into a 'prefix' for the melt
  line="";
  for (i=1; i<=9; i++)
    line=(line=="" ? $i : line OFS $i);

  # construct columns 12-15 into a 'suffix' for the melt
  suffline="";
  for (i=12; i<=15; i++)
    suffline=(suffline=="" ? $i : suffline OFS $i);

  # split exonstarts, exonends, and exonframes into array
  # this assumes they are in the same order
  StartCount = split($10,exStart,",");
  EndCount = split($11,exEnd,",");
  FrameCount = split($16,exFrame,",");

  if (StartCount-1 != $9 && exStart[$9] != "") {
    print "ExonCounts don't match with exonStarts" > "/dev/stderr"
    print "ExonCounts is: "$9" but found "StartCount" exons in line "NR > "/dev/stderr"
    print exStart[$9] > "/dev/stderr";
    print $0 > "/dev/stderr";
    exit 1;
  }
  
  if (EndCount-1 != $9 && exEnd[$9] != "")  {
    print "ExonCounts don't match with exonEnds" > "/dev/stderr";
    print "ExonCounts is: "$9" but found "EndCount" exons in line "NR > "/dev/stderr"
    print $0 > "/dev/stderr";
    exit 1;
  }

  if (FrameCount-1 != $9 && exFrame[$9] != "")  {
    print "ExonCounts don't match with exonFrames" > "/dev/stderr";
    print "ExonCounts is: "$9" but found "FrameCount" exons in line "NR > "/dev/stderr"
    print $0 > "/dev/stderr";
    exit 1;
  }
     
  # construct output 
  for (i=1; i<=$9; i++)
    print line OFS exStart[i] OFS exEnd[i] OFS suffline OFS exFrame[i];
}
