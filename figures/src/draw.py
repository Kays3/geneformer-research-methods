#!/usr/bin/env python3
"""Generate the methods diagrams offline: python figures/src/draw.py."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','font.size':11,'savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
INK='#17324D'; BLUE='#0072B2'; ORANGE='#D55E00'; GREEN='#00795E'; GREY='#586978'; PALE='#EDF4F8'; LIGHT='#FFF2E6'

def canvas(title, subtitle):
 fig,ax=plt.subplots(figsize=(12,6.4));fig.patch.set_facecolor('white');ax.set_xlim(0,12);ax.set_ylim(0,6.4);ax.axis('off')
 ax.text(.35,6.05,title,fontsize=21,fontweight='bold',color=INK,va='top');ax.text(.36,5.57,subtitle,fontsize=11,color=GREY,va='top');return fig,ax

def box(ax,x,y,w,h,title,body='',color=BLUE,fill=PALE):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.05,rounding_size=0.12',linewidth=1,edgecolor='#C5D3DB',facecolor=fill))
 ax.plot([x+.12,x+.12],[y+.14,y+h-.14],color=color,lw=3,solid_capstyle='round')
 ax.text(x+.27,y+h-.22,title,color=color,fontsize=11,fontweight='bold',va='top')
 if body:ax.text(x+.27,y+h-.57,body,color=INK,fontsize=10,va='top',linespacing=1.42)

def arrow(ax,x1,y1,x2,y2,color=GREY):
 ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=13,linewidth=1.4,color=color,connectionstyle='arc3'))

def save(fig,name):
 if len(sys.argv) > 1 and name not in sys.argv[1:]:
  plt.close(fig)
  return
 for ext in ('svg','png'):
  metadata = {'Creator': None, 'Date': None} if ext == 'svg' else {'Software': 'Matplotlib'}
  fig.savefig(OUT/f'{name}.{ext}',dpi=180,bbox_inches='tight',pad_inches=.12,metadata=metadata)
 plt.close(fig)

fig,ax=canvas('From counts to a reviewed readout','Each transition has an output and a check; model shifts are predictions.')
nodes=[('Prepare','Counts and\nprovenance'),('Tokenise','Ranked gene\nsequences'),('Embed','Released or\nfine-tuned\nmodel'),('Perturb','Delete or move\na token'),('Readout','Magnitude or\nstate direction'),('Test','Matched null;\ntest by donor')]
for i,(t,b) in enumerate(nodes):
 x=.35+i*1.92;box(ax,x,2.65,1.65,1.65,t,b,color=BLUE if i<3 else GREEN)
 if i<5:arrow(ax,x+1.72,3.47,x+1.88,3.47)
box(ax,.35,.8,5.45,1.05,'BEFORE','Fix the cohort, controls, rules and stop conditions.',ORANGE,LIGHT)
box(ax,6.15,.8,5.45,1.05,'AFTER','Record limits; science, numbers\nand privacy reviews.',GREEN)
save(fig,'01-pipeline')

fig,ax=canvas('An edit moves a representation along an axis','Illustrative geometry; points and arrow lengths are schematic.')
ax.annotate('',xy=(8.8,3.1),xytext=(2,3.1),arrowprops={'arrowstyle':'->','lw':2,'color':GREY});ax.text(1.8,2.72,'start centroid',color=INK);ax.text(8,2.72,'goal centroid',color=INK)
ax.scatter([3.9,4.85,4.35],[3.1,4.05,2.25],s=[110,100,90],c=[BLUE,ORANGE,GREEN],zorder=3)
arrow(ax,3.9,3.1,4.85,4.05,ORANGE);ax.plot([4.85,4.85],[3.1,4.05],ls=':',color=GREY);ax.annotate('',xy=(4.85,3.1),xytext=(3.9,3.1),arrowprops={'arrowstyle':'->','color':ORANGE,'lw':3})
ax.text(3.4,3.35,'original',color=BLUE);ax.text(4.98,4.02,'edited',color=ORANGE);ax.text(3.78,2.42,'matched edit',color=GREEN)
box(ax,8.55,3.95,3,1.25,'GOAL AXIS','Cross-fit and check\nresolution first.',BLUE)
box(ax,8.55,1.1,3,1.25,'MATCHED NULL','Compare with edits matched\non rank and detection.',GREEN)
ax.text(.4,.45,'An axis projection is not a measured effect in a cell.',color=GREY,fontsize=10)
save(fig,'02-axis-null')

fig,ax=canvas('Compare released and fine-tuned models on the same readouts','Keep donor groups disjoint and hold readout cells out of training.')
box(ax,.4,3.55,2.5,1.35,'DONOR SPLIT','Training · evaluation\n· test / held-out cells')
box(ax,3.65,4.1,3.1,1.05,'RELEASED MODEL','No model training.',BLUE)
box(ax,3.65,2.4,3.1,1.25,'FINE-TUNED MODEL','Train only on training donors.\nChoose on evaluation donors.',ORANGE,LIGHT)
box(ax,7.45,4.1,3.9,1.05,'MATCHED READOUTS','Same eligible cells and operations.',BLUE)
box(ax,7.45,2.4,3.9,1.25,'GATES BEFORE ISP','Classifier and collapse checks.\nStop the arm if a required gate fails.',ORANGE,LIGHT)
for a,b in [((2.95,4.25),(3.6,4.62)),((2.95,4.0),(3.6,3.03)),((6.8,4.62),(7.4,4.62)),((6.8,3.03),(7.4,3.03))]:arrow(ax,*a,*b)
arrow(ax,9.35,3.7,9.35,4.05,GREEN)
ax.text(.42,.9,'Centroids and axes are built for each model; held-out readout cells do not train the model.',color=INK,fontsize=11)
save(fig,'03-comparison-design')

fig,ax=canvas('Collapse can return a successful exit code','Illustration of diagnostics; no observed data or fitted values are shown.')
box(ax,.45,3.35,3.25,1.6,'SYMPTOMS','Nearly constant predictions\nAUC near chance\nEmbeddings nearly identical',ORANGE,LIGHT)
box(ax,4.35,3.35,3.25,1.6,'CHECK EVERY EPOCH','Evaluation AUC\nPrediction spread\nLoss and gradient logs',BLUE)
box(ax,8.25,3.35,3.25,1.6,'CHECK EMBEDDINGS','Pairwise cosine on a seeded\nstratified sample\nVoid a collapsed arm',GREEN)
arrow(ax,3.75,4.1,4.3,4.1);arrow(ax,7.65,4.1,8.2,4.1)
box(ax,.45,1.05,11.05,1.25,'INTERPRETATION','Cause not established. Keep the best eligible checkpoint; report failed epochs.\nA collapsed model can still emit numbers: those readings are void.',ORANGE,LIGHT)
save(fig,'04-collapse')

fig,ax=canvas('Register first, review independently','A changed input or rule requires a new pin and the appropriate re-review.')
for x,t,b,c in [(.35,'REGISTER','Question, null,\nfamily, stop rules',BLUE),(3.25,'FREEZE','Code, inputs,\ncontrols and hashes',GREEN),(6.15,'RUN','Start checks;\nrecord outputs',BLUE),(9.05,'REVIEW','Science · numbers\n· privacy',ORANGE)]:box(ax,x,2.95,2.5,1.55,t,b,c,LIGHT if c==ORANGE else PALE)
for x in (2.9,5.8,8.7):arrow(ax,x,3.72,x+.3,3.72)
ax.text(.42,1.45,'A gate checks evidence for its own scope. A clean self-check is not an independent pass.',color=INK,fontsize=11)
save(fig,'05-review-flow')

fig,ax=canvas('A practical route through common failures','Start at the symptom; confirm the cause before changing a run.')
box(ax,.35,3.35,3.05,1.45,'NO OR INVALID OUTPUT','Check completion, NaNs\nand the no-op control.',BLUE)
box(ax,4.45,3.90,3.0,1.25,'ALL SHIFTS ZERO?','Check precision and\nresolution.',ORANGE,LIGHT)
box(ax,4.45,2.35,3.0,1.1,'CONSTANT SCORES?','Check collapse diagnostics.',ORANGE,LIGHT)
box(ax,8.45,3.35,3.05,1.45,'REPAIR OR STOP','Trace the observed cause;\nkeep failed outputs and\nchange one thing per rerun.',GREEN)
for a,b in [((3.45,4.15),(4.4,4.55)),((3.45,3.95),(4.4,2.9)),((7.5,4.55),(8.4,4.15)),((7.5,2.9),(8.4,3.9))]:arrow(ax,*a,*b)
ax.text(.4,1.15,'A failed gate is not evidence of no biological effect.',color=GREY,fontsize=11)
save(fig,'06-troubleshooting')

rng=np.random.default_rng(173); controls=rng.normal(0,0.45,160); observations=np.array([-.32,.12,.26,.38,.55,.72,.13,.41])
fig,ax=plt.subplots(figsize=(12,6.4));fig.patch.set_facecolor('white');ax.set_facecolor('white')
fig.text(.12,.93,'A matched-null readout',fontsize=21,fontweight='bold',color=INK);fig.text(.12,.875,'SYNTHETIC EXAMPLE · invented values for a methods illustration',fontsize=11,color=ORANGE,fontweight='bold')
ax.hist(controls,bins=18,color='#C8DAE4',edgecolor='white',label='Matched edits (synthetic)');ax.axvline(np.median(observations),color=ORANGE,lw=2.5,label='Example median (synthetic)')
ax.set_xlabel('Illustrative shift (arbitrary units)',color=INK);ax.set_ylabel('Number of matched edits',color=INK);ax.legend(frameon=False,loc='upper left');ax.grid(axis='y',alpha=.16);fig.text(.12,.04,'The null is illustrative; this figure does not report a study result or p value.',color=GREY,fontsize=10)
fig.subplots_adjust(top=.78,bottom=.18,left=.12,right=.95);save(fig,'07-synthetic-readout')
