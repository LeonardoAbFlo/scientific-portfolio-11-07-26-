import glob, re, statistics as st

def avg_xvg(path):
  xs=[]
  with open(path) as f:
    for line in f:
      if line.startswith(('#','@')): 
        continue
      parts=line.split()
      if len(parts)>=2:
        xs.append(float(parts[1]))
  return sum(xs)/len(xs)

rmsd=[avg_xvg(p) for p in sorted(glob.glob("rep*/rmsd_rna.xvg"))]
rg=[avg_xvg(p) for p in sorted(glob.glob("rep*/rg_rna.xvg"))]
hb=[avg_xvg(p) for p in sorted(glob.glob("rep*/hb_stem.xvg"))]
ct=[avg_xvg(p) for p in sorted(glob.glob("rep*/contacts_stem.xvg"))]
ee=[avg_xvg(p) for p in sorted(glob.glob("rep*/end2end_com.xvg"))]

def pm_sd(arr):
  return (st.mean(arr), st.pstdev(arr))

for name, arr in [("RMSD (nm)",rmsd),("Rg (nm)",rg),("HB stem",hb),("Contacts stem",ct),("End2End (nm)",ee)]:
  m,s=pm_sd(arr)
  print(f"{name}: mean={m:.4f}, sd={s:.4f}, n={len(arr)}")
