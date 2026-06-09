from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
FIGS = ROOT / 'figures'
FIGS.mkdir(exist_ok=True)

# This script regenerates the traceability figures from the revision tables.
completion = pd.read_csv(DATA / 'reviewer_requirement_completion_audit_20260609.csv')
trace = pd.read_csv(DATA / 'traceability_total_evidence_matrix_20260609.csv')

concerns = ["c'/phi'\ndegradation", "eta(t)\nretention", "treated-zone\ncredit", "chi_s\nlocalization", "Fig. 2\ntime trend", "Fig. 4\ndesign map", "direct V5\nfooting"]
columns = ['Element\nstrength', 'Durability/\ndegradation', 'Footing/\nPLT', 'Field\nevidence', 'Numerical\naudit', 'Integrated\ndegraded footing']
Z = np.array([[3,2,1,0,1,0],[1,3,1,2,1,0],[1,0,3,2,3,0],[2,1,0,1,2,0],[2,3,1,1,2,0],[1,2,2,1,3,0],[0,0,0,0,0,0]], dtype=float)
plt.figure(figsize=(8.3,5.1))
im=plt.imshow(Z, aspect='auto')
plt.xticks(np.arange(len(columns)), columns, rotation=0, ha='center')
plt.yticks(np.arange(len(concerns)), concerns)
plt.title('Reviewer-requested evidence coverage map')
cbar=plt.colorbar(im)
cbar.set_ticks([0,1,2,3])
cbar.set_ticklabels(['absent', 'proxy', 'component', 'direct/near-direct'])
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        plt.text(j, i, ['0','P','C','D'][int(Z[i,j])], ha='center', va='center', fontsize=8)
plt.xlabel('Evidence class available in the resubmission package')
plt.tight_layout()
plt.savefig(FIGS / 'figure_traceability_reviewer_requirement_map.png', dpi=300)
plt.close()

components = ['c/phi parameters','degradation clock','treated-zone credit','localization chi_s','Figure 2 trend','Figure 4 map','direct V5 test']
status_labels = ['not identified','sensitivity','scenario-bounded','bounded','component-calibrated']
status_values = [4,3,3,1,2,2,0]
plt.figure(figsize=(8.2,4.7))
y = np.arange(len(components))
plt.scatter(status_values, y, s=90)
for val, yi in zip(status_values, y):
    plt.text(val+0.04, yi, status_labels[val], va='center', fontsize=8)
plt.yticks(y, components)
plt.xticks(range(len(status_labels)), status_labels, rotation=25, ha='right')
plt.xlim(-0.3,4.9)
plt.title('Claim status by model component (claim status, not predictive accuracy)')
plt.grid(axis='x', linestyle=':', linewidth=0.6)
plt.tight_layout()
plt.savefig(FIGS / 'figure_claim_status_by_component.png', dpi=300)
plt.close()

print('Traceability figures regenerated in', FIGS)
