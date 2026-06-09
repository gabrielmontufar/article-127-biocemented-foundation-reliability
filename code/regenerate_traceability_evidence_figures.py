
"""Regenerate traceability figures added for the GTENG-15832 AE response.
The script uses only CSV files distributed in the supplementary archive.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)
trace_df = pd.read_csv(DATA / "parameter_dataset_traceability_matrix_20260609.csv")
short = ['c,phi\ndegradation','treated zone\ncontribution','localization\nfragility','Fig.2\ntime capacity','Fig.4\ndesign map','experimental\ncomparison','direct V5\ngap']
features = ['Direct footing/bearing data','Strength parameter data','Durability/degradation data','Localization evidence','Numerical/reproducible checks']
M = np.array([[1,2,2,0,1],[2,0,1,0,2],[0,0,0,1,2],[1,0,2,0,2],[1,0,1,0,2],[2,1,2,0,2],[0,0,0,0,1]])
fig, ax = plt.subplots(figsize=(10.8,5.6), dpi=300)
im = ax.imshow(M, aspect='auto')
ax.set_xticks(range(len(features)), labels=features, rotation=35, ha='right')
ax.set_yticks(range(len(short)), labels=short)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, {0:'--',1:'proxy',2:'data'}[int(M[i,j])], ha='center', va='center', fontsize=8)
ax.set_title('Traceability matrix: reviewer concern -> dataset class -> claim boundary')
cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cbar.set_ticks([0,1,2]); cbar.set_ticklabels(['none','proxy/bound','data'])
fig.tight_layout(); fig.savefig(FIGS/'figure_traceability_evidence_map.png', dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(9.5,5.4), dpi=300)
y = np.arange(len(short)); ax.barh(y, trace_df['evidence_score_0_to_4'])
ax.set_yticks(y, labels=short); ax.set_xlim(0,4)
ax.set_xlabel('Evidence status score (0 = gap; 1 = sensitivity; 2 = bounded; 3 = component calibrated; 4 = direct integrated)')
ax.set_title('Reviewer-request completion status by model component')
for yi, score, stat in zip(y, trace_df['evidence_score_0_to_4'], trace_df['status']):
    ax.text(score + 0.05, yi, stat, va='center', fontsize=8)
fig.tight_layout(); fig.savefig(FIGS/'figure_reviewer_completion_status.png', dpi=300); plt.close(fig)

kul = pd.read_csv(DATA/'kulkarni2021_table6_bcr_srf_summary.csv')
fig, ax = plt.subplots(figsize=(9.5,5.4), dpi=300)
labels=[f"{g}\n{s}" for g,s in zip(kul['geometry'], kul['plate_size'])]
x=np.arange(len(kul)); ax.bar(x, kul['BCR'])
ax.set_xticks(x, labels=labels, rotation=35, ha='right')
ax.set_ylabel('Bearing capacity ratio, BCR')
ax.set_title('Direct short-term MICP plate-load evidence used to bound initial treated-zone contribution')
for xi, val in zip(x, kul['BCR']): ax.text(xi, val + 0.05, f'{val:.2f}', ha='center', va='bottom', fontsize=8)
fig.tight_layout(); fig.savefig(FIGS/'figure_kulkarni_bcr_traceability.png', dpi=300); plt.close(fig)

sharma = pd.read_csv(DATA/'sharma2021_digitized_durability_retention_data.csv')
fig, ax = plt.subplots(figsize=(8.5,5.2), dpi=300)
ax.plot(sharma['exposure_value'], sharma['retention'], marker='o', label='digitized/report retention')
if 'predicted_retention' in sharma.columns: ax.plot(sharma['exposure_value'], sharma['predicted_retention'], marker='s', label='fitted retention model')
for _,r in sharma.iterrows(): ax.annotate(r['role'], (r['exposure_value'], r['retention']), textcoords='offset points', xytext=(0,7), ha='center', fontsize=7)
ax.set_xlabel('Freeze-thaw exposure cycles'); ax.set_ylabel('Normalized UCS retention'); ax.set_ylim(0,1.05)
ax.set_title('Durability evidence used to bound degradation-retention scenarios'); ax.legend()
fig.tight_layout(); fig.savefig(FIGS/'figure_degradation_retention_traceability.png', dpi=300); plt.close(fig)

rows = ['c/phi degradation','treated-zone function','chi_s localization','time-capacity scenarios','depth-degradation map','direct V5 footing']
cols = ['data support','calibration/comparison','directness','universal claim']
G = np.array([[2,2,1,0],[2,2,1,0],[1,1,0,0],[2,2,1,0],[2,1,1,0],[0,0,0,0]])
fig, ax = plt.subplots(figsize=(9,5.2), dpi=300)
im = ax.imshow(G, aspect='auto')
ax.set_xticks(range(len(cols)), labels=cols, rotation=25, ha='right'); ax.set_yticks(range(len(rows)), labels=rows)
for i in range(G.shape[0]):
    for j in range(G.shape[1]): ax.text(j, i, {0:'gap',1:'bounded',2:'supported'}[int(G[i,j])], ha='center', va='center', fontsize=8)
ax.set_title('Closure map: what the revision substantiates versus what remains outside the claim')
cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cbar.set_ticks([0,1,2]); cbar.set_ticklabels(['gap','bounded','supported'])
fig.tight_layout(); fig.savefig(FIGS/'figure_validation_gap_closure_map.png', dpi=300); plt.close(fig)
print('Traceability figures regenerated successfully.')
