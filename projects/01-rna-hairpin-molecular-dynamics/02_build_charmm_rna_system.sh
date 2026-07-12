awk '
BEGIN{
  OFS=""
  apos = sprintf("%c",39)

  resmap["A"]=" RA"; resmap["C"]=" RC"; resmap["G"]=" RG"; resmap["U"]=" RU"
  atommap["OP1"]=" O1P"; atommap["OP2"]=" O2P"; atommap["OP3"]=" O3P"
}

$0 ~ /^(ATOM  |HETATM)/{
  atom = substr($0,13,4)
  gsub(/\*/, apos, atom)
  key  = atom
  gsub(/ /,"",key)
  if (key in atommap) atom = atommap[key]
  $0 = substr($0,1,12) atom substr($0,17)

  res  = substr($0,18,3)
  rkey = res
  gsub(/ /,"",rkey)
  if (rkey in resmap) res = resmap[rkey]
  $0 = substr($0,1,17) res substr($0,21)
}
{ print }
' RNA.pdb > RNA_charmm_fix.pdb
